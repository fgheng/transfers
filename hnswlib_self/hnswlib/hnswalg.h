#pragma once

#include "visited_list_pool.h"
#include "hnswlib.h"
#include <atomic>
#include <random>
#include <stdlib.h>
#include <assert.h>
#include <unordered_set>
#include <list>
#include <memory>
#include "./space_pq.h"

namespace hnswlib {
typedef unsigned int tableint;
typedef unsigned int linklistsizeint;

template<typename dist_t>
class HierarchicalNSW : public AlgorithmInterface<dist_t> {
 public:
    static const tableint MAX_LABEL_OPERATION_LOCKS = 65536;
    static const unsigned char DELETE_MARK = 0x01;

    size_t max_elements_{0};
    mutable std::atomic<size_t> cur_element_count{0};  // current number of elements
    size_t size_data_per_element_{0};
    size_t size_links_per_element_{0};
    mutable std::atomic<size_t> num_deleted_{0};  // number of deleted elements
    size_t M_{0};
    size_t maxM_{0};
    size_t maxM0_{0};
    size_t ef_construction_{0};
    size_t ef_{ 0 };

    double mult_{0.0}, revSize_{0.0};
    int maxlevel_{0};

    std::unique_ptr<VisitedListPool> visited_list_pool_{nullptr};

    // Locks operations with element by label value
    mutable std::vector<std::mutex> label_op_locks_;

    std::mutex global;
    std::vector<std::mutex> link_list_locks_;

    tableint enterpoint_node_{0};

    size_t size_links_level0_{0};
    size_t offsetData_{0}, offsetLevel0_{0}, label_offset_{ 0 };

    size_t data_level0_memory_size_{0};
    char *data_level0_memory_{nullptr};
    char **linkLists_{nullptr};
    std::vector<int> element_levels_;  // keeps level of each element

    size_t data_size_{0};

    DISTFUNC<dist_t> fstdistfunc_;
    void *dist_func_param_{nullptr};

    mutable std::mutex label_lookup_lock;  // lock for label_lookup_
    std::unordered_map<labeltype, tableint> label_lookup_;

    std::default_random_engine level_generator_;
    std::default_random_engine update_probability_generator_;

    mutable std::atomic<long> metric_distance_computations{0};
    mutable std::atomic<long> metric_hops{0};

    bool allow_replace_deleted_ = false;  // flag to replace deleted elements (marked as deleted) during insertions

    std::mutex deleted_elements_lock;  // lock for deleted_elements
    std::unordered_set<tableint> deleted_elements;  // contains internal ids of deleted elements

    int pq_M_ = 0;
    int pq_nbits_ = 0;
    int pq_dsub = 0;
    int pq_ks_ = 0;
    float scale_ = 1.0;
    float scale2_ = 1.0;

    std::vector<std::vector<float>> pq_centroids_;
    std::vector<float> pq_residuals_;


    HierarchicalNSW(SpaceInterface<dist_t> *s) {
    }


    HierarchicalNSW(
        SpaceInterface<dist_t> *s,
        const std::string &location,
        bool nmslib = false,
        size_t max_elements = 0,
        bool allow_replace_deleted = false)
        : allow_replace_deleted_(allow_replace_deleted) {
        loadIndex(location, s, max_elements);
    }


    HierarchicalNSW(
        SpaceInterface<dist_t> *s,
        size_t max_elements,
        size_t M = 16,
        size_t ef_construction = 200,
        size_t random_seed = 100,
        bool allow_replace_deleted = false)
        : label_op_locks_(MAX_LABEL_OPERATION_LOCKS),
            link_list_locks_(max_elements),
            element_levels_(max_elements),
            allow_replace_deleted_(allow_replace_deleted) {
        max_elements_ = max_elements;
        num_deleted_ = 0;
        data_size_ = s->get_data_size();
        fstdistfunc_ = s->get_dist_func();
        dist_func_param_ = s->get_dist_func_param();
        if ( M <= 10000 ) {
            M_ = M;
        } else {
            HNSWERR << "warning: M parameter exceeds 10000 which may lead to adverse effects." << std::endl;
            HNSWERR << "         Cap to 10000 will be applied for the rest of the processing." << std::endl;
            M_ = 10000;
        }
        maxM_ = M_;
        maxM0_ = M_ * 2;
        ef_construction_ = std::max(ef_construction, M_);
        ef_ = 10;

        level_generator_.seed(random_seed);
        update_probability_generator_.seed(random_seed + 1);

        size_links_level0_ = maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);
        size_data_per_element_ = size_links_level0_ + data_size_ + sizeof(labeltype);
        offsetData_ = size_links_level0_;
        label_offset_ = size_links_level0_ + data_size_;
        offsetLevel0_ = 0;

        data_level0_memory_size_ = max_elements_ * size_data_per_element_;
        data_level0_memory_ = (char *) malloc(data_level0_memory_size_);
        if (data_level0_memory_ == nullptr)
            throw std::runtime_error("Not enough memory");

        cur_element_count = 0;

        visited_list_pool_ = std::unique_ptr<VisitedListPool>(new VisitedListPool(1, max_elements));

        // initializations for special treatment of the first node
        enterpoint_node_ = -1;
        maxlevel_ = -1;

        linkLists_ = (char **) malloc(sizeof(void *) * max_elements_);
        if (linkLists_ == nullptr)
            throw std::runtime_error("Not enough memory: HierarchicalNSW failed to allocate linklists");
        size_links_per_element_ = maxM_ * sizeof(tableint) + sizeof(linklistsizeint);
        // mult_ = 1 / log(1.0 * M_);
        mult_ = 1 / log(1.0 * 4);
        revSize_ = 1.0 / mult_;
    }


    ~HierarchicalNSW() {
        clear();
    }

    void clear() {
        free(data_level0_memory_);
        data_level0_memory_ = nullptr;
        for (tableint i = 0; i < cur_element_count; i++) {
            if (element_levels_[i] > 0)
                free(linkLists_[i]);
        }
        free(linkLists_);
        linkLists_ = nullptr;
        cur_element_count = 0;
        visited_list_pool_.reset(nullptr);
    }


    struct CompareByFirst {
        constexpr bool operator()(std::pair<dist_t, tableint> const& a,
            std::pair<dist_t, tableint> const& b) const noexcept {
            return a.first < b.first;
        }
    };


    void setEf(size_t ef) {
        ef_ = ef;
    }


    inline std::mutex& getLabelOpMutex(labeltype label) const {
        // calculate hash
        size_t lock_id = label & (MAX_LABEL_OPERATION_LOCKS - 1);
        return label_op_locks_[lock_id];
    }


    inline labeltype getExternalLabel(tableint internal_id) const {
        labeltype return_label;
        memcpy(&return_label, (data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_), sizeof(labeltype));
        return return_label;
    }


    inline void setExternalLabel(tableint internal_id, labeltype label) const {
        memcpy((data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_), &label, sizeof(labeltype));
    }


    inline labeltype *getExternalLabeLp(tableint internal_id) const {
        return (labeltype *) (data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_);
    }


    size_t get_data_level0_memory_size() {
        return data_level0_memory_size_;
    }

    char* get_data_level0_memory() {
        return data_level0_memory_;
    }


    inline char *getDataByInternalId(tableint internal_id) const {
        return (data_level0_memory_ + internal_id * size_data_per_element_ + offsetData_);
    }


    int getRandomLevel(double reverse_size) {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double r = -log(distribution(level_generator_)) * reverse_size;
        return (int) r;
    }

    size_t getMaxElements() {
        return max_elements_;
    }

    size_t getCurrentElementCount() {
        return cur_element_count;
    }

    size_t getDeletedCount() {
        return num_deleted_;
    }

    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst>
    searchBaseLayer(tableint ep_id, const void *data_point, int layer) {
        VisitedList *vl = visited_list_pool_->getFreeVisitedList();
        vl_type *visited_array = vl->mass;
        vl_type visited_array_tag = vl->curV;

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidateSet;

        dist_t lowerBound;
        if (!isMarkedDeleted(ep_id)) {
            dist_t dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_, scale2_);
            top_candidates.emplace(dist, ep_id);
            lowerBound = dist;
            candidateSet.emplace(-dist, ep_id);
        } else {
            lowerBound = std::numeric_limits<dist_t>::max();
            candidateSet.emplace(-lowerBound, ep_id);
        }
        visited_array[ep_id] = visited_array_tag;

        while (!candidateSet.empty()) {
            std::pair<dist_t, tableint> curr_el_pair = candidateSet.top();
            if ((-curr_el_pair.first) > lowerBound && top_candidates.size() == ef_construction_) {
                break;
            }
            candidateSet.pop();

            tableint curNodeNum = curr_el_pair.second;

            std::unique_lock <std::mutex> lock(link_list_locks_[curNodeNum]);

            int *data;  // = (int *)(linkList0_ + curNodeNum * size_links_per_element0_);
            if (layer == 0) {
                data = (int*)get_linklist0(curNodeNum);
            } else {
                data = (int*)get_linklist(curNodeNum, layer);
//                    data = (int *) (linkLists_[curNodeNum] + (layer - 1) * size_links_per_element_);
            }
            size_t size = getListCount((linklistsizeint*)data);
            tableint *datal = (tableint *) (data + 1);
#ifdef USE_SSE
            _mm_prefetch((char *) (visited_array + *(data + 1)), _MM_HINT_T0);
            _mm_prefetch((char *) (visited_array + *(data + 1) + 64), _MM_HINT_T0);
            _mm_prefetch(getDataByInternalId(*datal), _MM_HINT_T0);
            _mm_prefetch(getDataByInternalId(*(datal + 1)), _MM_HINT_T0);
#endif

            for (size_t j = 0; j < size; j++) {
                tableint candidate_id = *(datal + j);
//                    if (candidate_id == 0) continue;
#ifdef USE_SSE
                _mm_prefetch((char *) (visited_array + *(datal + j + 1)), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*(datal + j + 1)), _MM_HINT_T0);
#endif
                if (visited_array[candidate_id] == visited_array_tag) continue;
                visited_array[candidate_id] = visited_array_tag;
                char *currObj1 = (getDataByInternalId(candidate_id));

                dist_t dist1 = fstdistfunc_(data_point, currObj1, dist_func_param_, scale2_);
                if (top_candidates.size() < ef_construction_ || lowerBound > dist1) {
                    candidateSet.emplace(-dist1, candidate_id);
#ifdef USE_SSE
                    _mm_prefetch(getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif

                    if (!isMarkedDeleted(candidate_id))
                        top_candidates.emplace(dist1, candidate_id);

                    if (top_candidates.size() > ef_construction_)
                        top_candidates.pop();

                    if (!top_candidates.empty())
                        lowerBound = top_candidates.top().first;
                }
            }
        }
        visited_list_pool_->releaseVisitedList(vl);

        return top_candidates;
    }


    // bare_bone_search means there is no check for deletions and stop condition is ignored in return of extra performance
    template <bool bare_bone_search = true, bool collect_metrics = false>
    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst>
    searchBaseLayerST(
        tableint ep_id,
        const void *data_point,
        size_t ef,
        float q_residual,
        BaseFilterFunctor* isIdAllowed = nullptr,
        BaseSearchStopCondition<dist_t>* stop_condition = nullptr) const {
        VisitedList *vl = visited_list_pool_->getFreeVisitedList();
        vl_type *visited_array = vl->mass;
        vl_type visited_array_tag = vl->curV;

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidate_set;

        dist_t lowerBound;
        if (bare_bone_search || 
            (!isMarkedDeleted(ep_id) && ((!isIdAllowed) || (*isIdAllowed)(getExternalLabel(ep_id))))) {
            char* ep_data = getDataByInternalId(ep_id);
            dist_t dist = fstdistfunc_(data_point, ep_data, dist_func_param_, scale2_);
                    // add residuals
                    // dist += q_residual;
                    // dist += pq_residuals_[*getExternalLabeLp(ep_id)];
            lowerBound = dist;
            top_candidates.emplace(dist, ep_id);
            if (!bare_bone_search && stop_condition) {
                stop_condition->add_point_to_result(getExternalLabel(ep_id), ep_data, dist);
            }
            candidate_set.emplace(-dist, ep_id);
        } else {
            lowerBound = std::numeric_limits<dist_t>::max();
            candidate_set.emplace(-lowerBound, ep_id);
        }

        visited_array[ep_id] = visited_array_tag;

        while (!candidate_set.empty()) {
            std::pair<dist_t, tableint> current_node_pair = candidate_set.top();
            dist_t candidate_dist = -current_node_pair.first;

            bool flag_stop_search;
            if (bare_bone_search) {
                flag_stop_search = candidate_dist > lowerBound;
            } else {
                if (stop_condition) {
                    flag_stop_search = stop_condition->should_stop_search(candidate_dist, lowerBound);
                } else {
                    flag_stop_search = candidate_dist > lowerBound && top_candidates.size() == ef;
                }
            }
            if (flag_stop_search) {
                break;
            }
            candidate_set.pop();

            tableint current_node_id = current_node_pair.second;
            int *data = (int *) get_linklist0(current_node_id);
            size_t size = getListCount((linklistsizeint*)data);
//                bool cur_node_deleted = isMarkedDeleted(current_node_id);
            if (collect_metrics) {
                metric_hops++;
                metric_distance_computations+=size;
            }

#ifdef USE_SSE
            _mm_prefetch((char *) (visited_array + *(data + 1)), _MM_HINT_T0);
            _mm_prefetch((char *) (visited_array + *(data + 1) + 64), _MM_HINT_T0);
            _mm_prefetch(data_level0_memory_ + (*(data + 1)) * size_data_per_element_ + offsetData_, _MM_HINT_T0);
            _mm_prefetch((char *) (data + 2), _MM_HINT_T0);
#endif

            for (size_t j = 1; j <= size; j++) {
                int candidate_id = *(data + j);
//                    if (candidate_id == 0) continue;
#ifdef USE_SSE
                _mm_prefetch((char *) (visited_array + *(data + j + 1)), _MM_HINT_T0);
                _mm_prefetch(data_level0_memory_ + (*(data + j + 1)) * size_data_per_element_ + offsetData_,
                                _MM_HINT_T0);  ////////////
#endif
                if (!(visited_array[candidate_id] == visited_array_tag)) {
                    visited_array[candidate_id] = visited_array_tag;

                    char *currObj1 = (getDataByInternalId(candidate_id));
                    dist_t dist = fstdistfunc_(data_point, currObj1, dist_func_param_, scale2_);
                    // add residuals
                    // dist += q_residual;
                    // dist += pq_residuals_[*getExternalLabeLp(candidate_id)];

                    bool flag_consider_candidate;
                    if (!bare_bone_search && stop_condition) {
                        flag_consider_candidate = stop_condition->should_consider_candidate(dist, lowerBound);
                    } else {
                        flag_consider_candidate = top_candidates.size() < ef || lowerBound > dist;
                    }

                    if (flag_consider_candidate) {
                        candidate_set.emplace(-dist, candidate_id);
#ifdef USE_SSE
                        _mm_prefetch(data_level0_memory_ + candidate_set.top().second * size_data_per_element_ +
                                        offsetLevel0_,  ///////////
                                        _MM_HINT_T0);  ////////////////////////
#endif

                        if (bare_bone_search || 
                            (!isMarkedDeleted(candidate_id) && ((!isIdAllowed) || (*isIdAllowed)(getExternalLabel(candidate_id))))) {
                            top_candidates.emplace(dist, candidate_id);
                            if (!bare_bone_search && stop_condition) {
                                stop_condition->add_point_to_result(getExternalLabel(candidate_id), currObj1, dist);
                            }
                        }

                        bool flag_remove_extra = false;
                        if (!bare_bone_search && stop_condition) {
                            flag_remove_extra = stop_condition->should_remove_extra();
                        } else {
                            flag_remove_extra = top_candidates.size() > ef;
                        }
                        while (flag_remove_extra) {
                            tableint id = top_candidates.top().second;
                            top_candidates.pop();
                            if (!bare_bone_search && stop_condition) {
                                stop_condition->remove_point_from_result(getExternalLabel(id), getDataByInternalId(id), dist);
                                flag_remove_extra = stop_condition->should_remove_extra();
                            } else {
                                flag_remove_extra = top_candidates.size() > ef;
                            }
                        }

                        if (!top_candidates.empty())
                            lowerBound = top_candidates.top().first;
                    }
                }
            }
        }

        visited_list_pool_->releaseVisitedList(vl);
        return top_candidates;
    }


    void getNeighborsByHeuristic2(
        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> &top_candidates,
        const size_t M) {

        if (top_candidates.size() < M) {
            return;
        }

        std::priority_queue<std::pair<dist_t, tableint>> queue_closest;
        std::vector<std::pair<dist_t, tableint>> return_list;
        while (top_candidates.size() > 0) {
            queue_closest.emplace(-top_candidates.top().first, top_candidates.top().second);
            top_candidates.pop();
        }

        while (queue_closest.size()) {
            if (return_list.size() >= M)
                break;
            std::pair<dist_t, tableint> curent_pair = queue_closest.top();
            dist_t dist_to_query = -curent_pair.first;
            queue_closest.pop();
            bool good = true;

            for (std::pair<dist_t, tableint> second_pair : return_list) {
                dist_t curdist =
                        fstdistfunc_(getDataByInternalId(second_pair.second),
                                        getDataByInternalId(curent_pair.second),
                                        dist_func_param_, scale2_);
                if (curdist < dist_to_query) {
                    good = false;
                    break;
                }
            }
            if (good) {
                return_list.push_back(curent_pair);
            }
        }

        for (std::pair<dist_t, tableint> curent_pair : return_list) {
            top_candidates.emplace(-curent_pair.first, curent_pair.second);
        }
    }


    linklistsizeint *get_linklist0(tableint internal_id) const {
        return (linklistsizeint *) (data_level0_memory_ + internal_id * size_data_per_element_ + offsetLevel0_);
    }


    linklistsizeint *get_linklist0(tableint internal_id, char *data_level0_memory_) const {
        return (linklistsizeint *) (data_level0_memory_ + internal_id * size_data_per_element_ + offsetLevel0_);
    }


    linklistsizeint *get_linklist(tableint internal_id, int level) const {
        return (linklistsizeint *) (linkLists_[internal_id] + (level - 1) * size_links_per_element_);
    }


    linklistsizeint *get_linklist_at_level(tableint internal_id, int level) const {
        return level == 0 ? get_linklist0(internal_id) : get_linklist(internal_id, level);
    }


    tableint mutuallyConnectNewElement(
        const void *data_point,
        tableint cur_c,
        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> &top_candidates,
        int level,
        bool isUpdate,
        bool use_heuristic2 = true) {
        size_t Mcurmax = level ? maxM_ : maxM0_;
        if (use_heuristic2) getNeighborsByHeuristic2(top_candidates, M_);  // 启发式算法找到 M 个邻居
        if (use_heuristic2 && top_candidates.size() > M_)
            throw std::runtime_error("Should be not be more than M_ candidates returned by the heuristic");

        std::vector<tableint> selectedNeighbors;
        selectedNeighbors.reserve(M_);
        while (top_candidates.size() > 0) {
            selectedNeighbors.push_back(top_candidates.top().second);
            top_candidates.pop();
        }

        tableint next_closest_entry_point = selectedNeighbors.back();

        if (use_heuristic2)
            //
        {
            // lock only during the update
            // because during the addition the lock for cur_c is already acquired
            std::unique_lock <std::mutex> lock(link_list_locks_[cur_c], std::defer_lock);
            if (isUpdate) {
                lock.lock();
            }
            linklistsizeint *ll_cur;
            if (level == 0)
                ll_cur = get_linklist0(cur_c); // 获取 level0 的数据，包括邻居，数据，外部 id
            else
                ll_cur = get_linklist(cur_c, level);  // 获取 level 层的邻居数据

            if (*ll_cur && !isUpdate) {
                throw std::runtime_error("The newly inserted element should have blank link list");
            }
            setListCount(ll_cur, selectedNeighbors.size());
            tableint *data = (tableint *) (ll_cur + 1);
            for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {
                if (data[idx] && !isUpdate)
                    throw std::runtime_error("Possible memory corruption");
                if (level > element_levels_[selectedNeighbors[idx]]) // 判断当前点要插入的层是否超过了其邻居所在的最高层
                    throw std::runtime_error("Trying to make a link on a non-existent level");

                data[idx] = selectedNeighbors[idx]; // 存储邻居 
            }
        }

        for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {
            std::unique_lock <std::mutex> lock(link_list_locks_[selectedNeighbors[idx]]);

            linklistsizeint *ll_other;
            if (level == 0)
                ll_other = get_linklist0(selectedNeighbors[idx]);
            else
                ll_other = get_linklist(selectedNeighbors[idx], level); // 获取已经连接的邻居的数据 

            size_t sz_link_list_other = getListCount(ll_other);

            if (sz_link_list_other > Mcurmax)
                throw std::runtime_error("Bad value of sz_link_list_other");
            if (selectedNeighbors[idx] == cur_c)
                throw std::runtime_error("Trying to connect an element to itself");
            if (level > element_levels_[selectedNeighbors[idx]])  // 当前层是否超过了邻居所在的层，是的话邻居就不需要在当前点反向连接了
                                                                  // 只会在小于邻居的层进行互联，即使这个邻居没有在该层
                throw std::runtime_error("Trying to make a link on a non-existent level");

            tableint *data = (tableint *) (ll_other + 1); // 获取已经连接的邻居的邻居信息 

            bool is_cur_c_present = false;
            if (isUpdate) {
                for (size_t j = 0; j < sz_link_list_other; j++) {
                    if (data[j] == cur_c) {
                        is_cur_c_present = true;
                        break;
                    }
                }
            }

            if (use_heuristic2) {
            // If cur_c is already present in the neighboring connections of `selectedNeighbors[idx]` then no need to modify any connections or run the heuristics.
            if (!is_cur_c_present) {
                if (sz_link_list_other < Mcurmax) { // 如果邻居的邻居数量少于 M 个，那么直接将当前点插入到邻居的邻居列表中，构成双向图
                    data[sz_link_list_other] = cur_c;
                    setListCount(ll_other, sz_link_list_other + 1);
                } else {  
                    // 如果邻居的邻居数量已经超过了 M 个，那么使用启发式算法，重新从 M+1 个邻居中选择 M 个邻居，
                    // 存在当前点不是邻居最合适的点的请，所以也就导致了hnsw 图不一定是一个完全的双向图
                    // finding the "weakest" element to replace it with the new one
                    dist_t d_max = fstdistfunc_(getDataByInternalId(cur_c), getDataByInternalId(selectedNeighbors[idx]),
                                                dist_func_param_, scale2_);
                    // Heuristic:
                    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidates;
                    candidates.emplace(d_max, cur_c);

                    for (size_t j = 0; j < sz_link_list_other; j++) {
                        candidates.emplace(
                                fstdistfunc_(getDataByInternalId(data[j]), getDataByInternalId(selectedNeighbors[idx]),
                                                dist_func_param_, scale2_), data[j]);
                    }

                    getNeighborsByHeuristic2(candidates, Mcurmax);

                    int indx = 0;
                    while (candidates.size() > 0) {
                        data[indx] = candidates.top().second;
                        candidates.pop();
                        indx++;
                    }

                    setListCount(ll_other, indx);
                    // Nearest K:
                    /*int indx = -1;
                    for (int j = 0; j < sz_link_list_other; j++) {
                        dist_t d = fstdistfunc_(getDataByInternalId(data[j]), getDataByInternalId(rez[idx]), dist_func_param_);
                        if (d > d_max) {
                            indx = j;
                            d_max = d;
                        }
                    }
                    if (indx >= 0) {
                        data[indx] = cur_c;
                    } */
                }
            }

            } else {
                if (sz_link_list_other < Mcurmax) { // 如果邻居的邻居数量少于 M 个，那么直接将当前点插入到邻居的邻居列表中，构成双向图
                    data[sz_link_list_other] = cur_c;
                    setListCount(ll_other, sz_link_list_other + 1);
                }

            }

        }

        return next_closest_entry_point;
    }


    void resizeIndex(size_t new_max_elements) {
        if (new_max_elements < cur_element_count)
            throw std::runtime_error("Cannot resize, max element is less than the current number of elements");

        visited_list_pool_.reset(new VisitedListPool(1, new_max_elements));

        element_levels_.resize(new_max_elements);

        std::vector<std::mutex>(new_max_elements).swap(link_list_locks_);

        // Reallocate base layer
        char * data_level0_memory_new = (char *) realloc(data_level0_memory_, new_max_elements * size_data_per_element_);
        if (data_level0_memory_new == nullptr)
            throw std::runtime_error("Not enough memory: resizeIndex failed to allocate base layer");
        data_level0_memory_ = data_level0_memory_new;

        // Reallocate all other layers
        char ** linkLists_new = (char **) realloc(linkLists_, sizeof(void *) * new_max_elements);
        if (linkLists_new == nullptr)
            throw std::runtime_error("Not enough memory: resizeIndex failed to allocate other layers");
        linkLists_ = linkLists_new;

        max_elements_ = new_max_elements;
    }

    size_t indexFileSize() const {
        size_t size = 0;
        size += sizeof(offsetLevel0_);
        size += sizeof(max_elements_);
        size += sizeof(cur_element_count);
        size += sizeof(size_data_per_element_);
        size += sizeof(label_offset_);
        size += sizeof(offsetData_);
        size += sizeof(maxlevel_);
        size += sizeof(enterpoint_node_);
        size += sizeof(maxM_);

        size += sizeof(maxM0_);
        size += sizeof(M_);
        size += sizeof(mult_);
        size += sizeof(ef_construction_);

        size += cur_element_count * size_data_per_element_;

        for (size_t i = 0; i < cur_element_count; i++) {
            unsigned int linkListSize = element_levels_[i] > 0 ? size_links_per_element_ * element_levels_[i] : 0;
            size += sizeof(linkListSize);
            size += linkListSize;
        }
        return size;
    }

    void saveIndex(const std::string &location) {
        std::cout << "max level: " << maxlevel_ << std::endl;
        std::ofstream output(location, std::ios::binary);
        std::streampos position;

        writeBinaryPOD(output, offsetLevel0_);
        writeBinaryPOD(output, max_elements_);
        writeBinaryPOD(output, cur_element_count);
        writeBinaryPOD(output, size_data_per_element_);
        writeBinaryPOD(output, label_offset_);
        writeBinaryPOD(output, offsetData_);
        writeBinaryPOD(output, maxlevel_);
        writeBinaryPOD(output, enterpoint_node_);
        writeBinaryPOD(output, maxM_);

        writeBinaryPOD(output, maxM0_);
        writeBinaryPOD(output, M_);
        writeBinaryPOD(output, mult_);
        writeBinaryPOD(output, ef_construction_);

        output.write(data_level0_memory_, cur_element_count * size_data_per_element_);

        for (size_t i = 0; i < cur_element_count; i++) {
            unsigned int linkListSize = element_levels_[i] > 0 ? size_links_per_element_ * element_levels_[i] : 0;
            writeBinaryPOD(output, linkListSize);
            if (linkListSize)
                output.write(linkLists_[i], linkListSize);
        }
        output.close();
    }


    void loadIndex(const std::string &location, SpaceInterface<dist_t> *s, size_t max_elements_i = 0) {
        std::ifstream input(location, std::ios::binary);

        if (!input.is_open())
            throw std::runtime_error("Cannot open file");

        clear();
        // get file size:
        input.seekg(0, input.end);
        std::streampos total_filesize = input.tellg();
        input.seekg(0, input.beg);

        readBinaryPOD(input, offsetLevel0_);
        readBinaryPOD(input, max_elements_);
        readBinaryPOD(input, cur_element_count);

        size_t max_elements = max_elements_i;
        if (max_elements < cur_element_count)
            max_elements = max_elements_;
        max_elements_ = max_elements;
        readBinaryPOD(input, size_data_per_element_);
        readBinaryPOD(input, label_offset_);
        readBinaryPOD(input, offsetData_);
        readBinaryPOD(input, maxlevel_);
        readBinaryPOD(input, enterpoint_node_);
        std::cout << "++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "load index max level: " << maxlevel_ << std::endl;

        readBinaryPOD(input, maxM_);
        readBinaryPOD(input, maxM0_);
        readBinaryPOD(input, M_);
        readBinaryPOD(input, mult_);
        readBinaryPOD(input, ef_construction_);

        data_size_ = s->get_data_size();
        fstdistfunc_ = s->get_dist_func();
        dist_func_param_ = s->get_dist_func_param();

        auto pos = input.tellg();

        /// Optional - check if index is ok:
        input.seekg(cur_element_count * size_data_per_element_, input.cur);
        for (size_t i = 0; i < cur_element_count; i++) {
            if (input.tellg() < 0 || input.tellg() >= total_filesize) {
                throw std::runtime_error("Index seems to be corrupted or unsupported");
            }

            unsigned int linkListSize;
            readBinaryPOD(input, linkListSize);
            if (linkListSize != 0) {
                input.seekg(linkListSize, input.cur);
            }
        }

        // throw exception if it either corrupted or old index
        if (input.tellg() != total_filesize)
            throw std::runtime_error("Index seems to be corrupted or unsupported");

        input.clear();
        /// Optional check end

        input.seekg(pos, input.beg);

        data_level0_memory_ = (char *) malloc(max_elements * size_data_per_element_);
        if (data_level0_memory_ == nullptr)
            throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");
        input.read(data_level0_memory_, cur_element_count * size_data_per_element_);

        size_links_per_element_ = maxM_ * sizeof(tableint) + sizeof(linklistsizeint);

        size_links_level0_ = maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);
        std::vector<std::mutex>(max_elements).swap(link_list_locks_);
        std::vector<std::mutex>(MAX_LABEL_OPERATION_LOCKS).swap(label_op_locks_);

        visited_list_pool_.reset(new VisitedListPool(1, max_elements));

        linkLists_ = (char **) malloc(sizeof(void *) * max_elements);
        if (linkLists_ == nullptr)
            throw std::runtime_error("Not enough memory: loadIndex failed to allocate linklists");
        element_levels_ = std::vector<int>(max_elements);
        revSize_ = 1.0 / mult_;
        ef_ = 10;
        for (size_t i = 0; i < cur_element_count; i++) {
            label_lookup_[getExternalLabel(i)] = i;
            unsigned int linkListSize;
            readBinaryPOD(input, linkListSize);
            if (linkListSize == 0) {
                element_levels_[i] = 0;
                linkLists_[i] = nullptr;
            } else {
                element_levels_[i] = linkListSize / size_links_per_element_;
                linkLists_[i] = (char *) malloc(linkListSize);
                if (linkLists_[i] == nullptr)
                    throw std::runtime_error("Not enough memory: loadIndex failed to allocate linklist");
                input.read(linkLists_[i], linkListSize);
            }
        }

        for (size_t i = 0; i < cur_element_count; i++) {
            if (isMarkedDeleted(i)) {
                num_deleted_ += 1;
                if (allow_replace_deleted_) deleted_elements.insert(i);
            }
        }

        input.close();

        return;
    }

    int getElementLevel(tableint internal_id) {
        return element_levels_[internal_id];
    }

    void getExternalNeighbours(tableint internal_id, std::vector<labeltype>& external_neighbours) const {
        linklistsizeint *ll_cur = get_linklist0(internal_id);
        tableint *data = (tableint *) (ll_cur + 1);

        for (int i = 0; i < *ll_cur; i++) {
            tableint internal_id_neighbour = *(data+i);

            labeltype external_id = getExternalLabel(internal_id_neighbour);
            external_neighbours.push_back(external_id);
        }
    }

    void countOutDegrees(std::vector<std::vector<linklistsizeint>>& out_degrees) {
        out_degrees.resize(max_elements_);
        for (int i = 0; i < max_elements_; i++) {
            out_degrees[i].resize(element_levels_[i]+1);
        }

        for (int i = 0; i < max_elements_; i++) {
            for (int level = 0; level <= element_levels_[i]; level++) {
                linklistsizeint* ll_cur;
                if (level == 0) {
                    ll_cur = get_linklist0(i);
                } else {
                    ll_cur = get_linklist(i, level);
                }
                out_degrees[i][level] = *ll_cur;
            }
        }
    }

    void countInDegrees(std::vector<std::vector<linklistsizeint>>& in_degrees) {
        in_degrees.resize(max_elements_);
        for (int i = 0; i < max_elements_; i++) {
            in_degrees[i].resize(element_levels_[i]+1);
        }

        for (int i = 0; i < max_elements_; i++) {
            for (int level = 0; level <= element_levels_[i]; level++) {
                linklistsizeint* ll_cur;
                if (level == 0) {
                    ll_cur = get_linklist0(i);
                } else {
                    ll_cur = get_linklist(i, level);
                }

                tableint *data = (tableint *) (ll_cur + 1);
                for (int j = 0; j < *ll_cur; j++) {
                    tableint internal_id_neighbour = *(data+j);
                    in_degrees[internal_id_neighbour][level]++;
                }
            }
        }
    }

    // 索引合并测试
    // 仅测试召回功能，不涉及内存优化
    void mergeIndex(const std::vector<HierarchicalNSW*>& shard_indexes) {
        if (shard_indexes.size() <= 1) {
            return;
        }

        struct node {
            labeltype external_label_;
            tableint internal_label_;
            uint32_t shard_id_;
            int max_level_;
            std::vector<std::vector<labeltype>> external_neighbours_; // 记录每个 level 的邻居
            std::vector<std::vector<tableint>> internal_neighbours_;
        };

        // 拿到所有点每层的边关系，构造自定义图
        std::cout << "get all edges from shards" << std::endl;
        std::vector<struct node> graph;
        for (int shard_id = 0; shard_id < shard_indexes.size(); shard_id++) {
            auto& index = shard_indexes[shard_id];
            size_t shard_max_elements = index->max_elements_;

            for (int shard_internal_label = 0; shard_internal_label < shard_max_elements; shard_internal_label++) {
                int cur_level = index->getElementLevel(shard_internal_label);

                struct node node;
                node.max_level_ = cur_level;
                node.external_label_ = index->getExternalLabel(shard_internal_label);
                node.internal_label_ = shard_internal_label;
                node.shard_id_ = shard_id;

                for (int level = 0; level <= cur_level; level++) {
                    std::vector<labeltype> neighbours;
                    if (level == 0) {
                        index->getExternalNeighbours(shard_internal_label, neighbours);
                    } else {
                        linklistsizeint* ll_cur = index->get_linklist(shard_internal_label, level);
                        tableint* data = (tableint*)(ll_cur+1);

                        // 遍历当前层的邻居
                        for (int i = 0; i < *ll_cur; i++) {
                            tableint internal_id_neighbour = *(data+i);
                            labeltype external_id = index->getExternalLabel(internal_id_neighbour);
                            neighbours.push_back(external_id);
                        }
                    }
                    node.external_neighbours_.push_back(neighbours);
                }

                graph.push_back(node);
            }
        }

        // 遍历所有的 external_label，将数据复制到当前索引中
        // 主要存储向量和外部标签，实际的邻居数量以及实际的邻居需要合并后再填充
        if (graph.empty()) {
            throw std::runtime_error("No edges to merge");
        }

        // 重排，将一样的 external_label 的边放在一起 
        std::cout << "resort graph nodes" << std::endl;
        std::sort(graph.begin(), graph.end(), [](const auto &left, const auto &right) {
            return left.external_label_ < right.external_label_ || 
                    (left.external_label_ == right.external_label_ && left.shard_id_ < right.shard_id_);
        });

        std::cout << "cal unique_labels size" << std::endl;
        size_t total_unique_elements = 0;
        std::unordered_set<labeltype> unique_labels;
        for (const auto& node : graph) {
            unique_labels.insert(node.external_label_);
        }
        total_unique_elements = unique_labels.size();
        std::cout << "unique_labels size: " << unique_labels.size() << std::endl;
        std::cout << "graph nodes size: " << graph.size() << std::endl;

        cur_element_count = total_unique_elements;

        // 合并边，存储到自定义图中
        std::cout << "merge edges" << std::endl;
        std::vector<struct node> merge_graph(total_unique_elements);  // 存储合并后的图
        std::unordered_map<labeltype, tableint> external_label_to_internal_id;  // 记录填充后的外部数据在当前索引中的位置
        labeltype current_external_label = graph[0].external_label_;
        tableint current_internal_label = 0;
        int shard_id = graph[0].shard_id_;
        auto& shard_index = shard_indexes[shard_id];
        int cur_max_level = graph[0].max_level_;

        // 合并 0 号点的边
        // std::vector<std::vector<std::vector<labeltype>>> graph_merge_neighbours(max_elements_);
        std::vector<std::vector<labeltype>> merge_neighbours(cur_max_level+1);
        // 不用临时变量，直接将 merge_graph 里的external_label_大小更改，直接向 merge_graph 里添加数据
        std::cout << "resize merge_neighbours size: " << std::endl;
        cur_max_level = graph[0].max_level_;
        for (int i = 1; i < graph.size(); i++) {
            if (graph[i].external_label_ == current_external_label) {
                if (graph[i].max_level_ > cur_max_level) {
                    cur_max_level = graph[i].max_level_;
                }
            } else {
                merge_graph[current_internal_label].external_neighbours_.resize(cur_max_level+1);
                merge_graph[current_internal_label].internal_neighbours_.resize(cur_max_level+1);
                merge_graph[current_internal_label].max_level_ = cur_max_level;
                merge_graph[current_internal_label].internal_label_ = current_internal_label;
                merge_graph[current_internal_label].external_label_ = current_external_label;
                current_internal_label++;

                cur_max_level = graph[i].max_level_;
                current_external_label = graph[i].external_label_;
            }

            // 更新 maxlevel_
            if (cur_max_level > maxlevel_) {
                maxlevel_ = cur_max_level;
            }
        }
        std::cout << "max level: " << maxlevel_ << std::endl;
        merge_graph[current_internal_label].external_neighbours_.resize(cur_max_level+1);
        merge_graph[current_internal_label].internal_neighbours_.resize(cur_max_level+1);
        merge_graph[current_internal_label].max_level_ = cur_max_level;
        merge_graph[current_internal_label].internal_label_ = current_internal_label;
        merge_graph[current_internal_label].external_label_ = current_external_label;

        current_internal_label = 0;
        current_external_label = graph[0].external_label_;

        for (int level = 0; level <= graph[0].max_level_; level++) {
            auto& merge_graph_neighbours_level  = merge_graph[current_internal_label].external_neighbours_[level];
            auto& graph_external_neighbours_level = graph[0].external_neighbours_[level];
            merge_graph_neighbours_level.insert(merge_graph_neighbours_level.end(), graph_external_neighbours_level.begin(), graph_external_neighbours_level.end());
        }

        std::cout << "merge all edges" << std::endl;
        std::cout << "graph size: " << graph.size() << std::endl;

        memcpy(getDataByInternalId(current_internal_label), 
               shard_index->getDataByLabelFloat(current_external_label).data(), data_size_);  // 复制向量
        memcpy(getExternalLabeLp(current_internal_label),
               &current_external_label, sizeof(labeltype));                     // 复制外部id
        external_label_to_internal_id[current_external_label] = current_internal_label;    // 记录外部id在当前索引中的位置 

        // 遍历所有的边，合并各个 level 的邻居
        for (int i = 1; i < graph.size(); i++) {
            if (graph[i].external_label_ != current_external_label) {
                current_internal_label++;
                current_external_label = graph[i].external_label_;

                // 将数据 copy 到 level0
                uint32_t shard_id = graph[i].shard_id_;
                auto& shard_index = shard_indexes[shard_id];
                memcpy(getDataByInternalId(current_internal_label), 
                    // shard_index->template getDataByLabel<float>(current_external_label).data(), data_size_);  // 复制向量
                    shard_index->getDataByLabelFloat(current_external_label).data(), data_size_);  // 复制向量
                memcpy(getExternalLabeLp(current_internal_label),
                    &current_external_label, sizeof(labeltype)); // 复制外部id
                external_label_to_internal_id[current_external_label] = current_internal_label;  // 记录外部id在当前索引中的位置 
                                                                                                 //
            }

            for (int level = 0; level <= graph[i].max_level_; level++) {
                auto& merge_graph_neighbours_level  = merge_graph[current_internal_label].external_neighbours_[level];
                auto& graph_external_neighbours_level = graph[i].external_neighbours_[level];
                merge_graph_neighbours_level.insert(merge_graph_neighbours_level.end(), graph_external_neighbours_level.begin(), graph_external_neighbours_level.end());
            }
        }

        // 邻居去重
        std::cout << "unique" << std::endl;
        for (auto& node: merge_graph) {
            auto& external_neighbours = node.external_neighbours_;
            for (int level = 0; level <= node.max_level_; level++) {
                auto& merge_neighbours_level = external_neighbours[level];
                std::sort(merge_neighbours_level.begin(), merge_neighbours_level.end());
                auto it = std::unique(merge_neighbours_level.begin(), merge_neighbours_level.end());
                merge_neighbours_level.erase(it, merge_neighbours_level.end());
            }
        }

        // 检查 max_level跟实际的邻居数量是否一致 
        std::cout << "merge graph size: " << merge_graph.size() << std::endl;
        for (int i = 0; i < merge_graph.size(); i++) {
            if (merge_graph[i].max_level_ != merge_graph[i].external_neighbours_.size() - 1) {
                std::cout << "The edges are not consistent: max_level: " << merge_graph[i].max_level_ << " neighbours size: " << merge_graph[i].external_neighbours_.size() << std::endl;
                throw std::runtime_error("The edges are not consistent");
            }
        }

        // 添加反向边
        // std::cout << "reconnect" << std::endl;
        // for (int i = 0; i < merge_graph.size(); i++) {
        //     auto external_label = merge_graph[i].external_label_;
        //     auto& external_neighbours = merge_graph[i].external_neighbours_;
        //     for (int level = 0; level <= merge_graph[i].max_level_; level++) {
        //         auto& external_neighbours_level = external_neighbours[level];
        //         for (int j = 0; j < external_neighbours_level.size(); j++) {
        //             auto neighbour_external_id = external_neighbours_level[j];
        //             auto neighbour_internal_id = external_label_to_internal_id.find(neighbour_external_id);
        //
        //             auto& neighbour_node = merge_graph[neighbour_internal_id->second];
        //             auto& neighbour_node_external_neighbours_level = neighbour_node.external_neighbours_[level];
        //             bool update = true;
        //             for (int k = 0; k < neighbour_node_external_neighbours_level.size(); k++) {
        //                 if (neighbour_node_external_neighbours_level[k] == external_label) {
        //                     update = false;
        //                     break;
        //                 }
        //             }
        //
        //             if (update) {
        //                 neighbour_node_external_neighbours_level.push_back(external_label);
        //             }
        //         }
        //     }
        // }


        // 映射内部 id
        std::cout << "map internal id" << std::endl;
        for (auto& node: merge_graph) {
            auto& external_neighbours = node.external_neighbours_;
            auto& internal_neighbours = node.internal_neighbours_;

            for (int level = 0; level <= node.max_level_; level++) {
                internal_neighbours[level].resize(external_neighbours[level].size());
            }

            for (int level = 0; level <= node.max_level_; level++) {
                for (int i = 0; i < external_neighbours[level].size(); i++) {
                    auto it = external_label_to_internal_id.find(external_neighbours[level][i]);
                    if (it == external_label_to_internal_id.end()) {
                        std::cout << "error: external id not found" << std::endl;
                        throw std::runtime_error("external id not found");
                    } else {
                        internal_neighbours[level][i] = (it->second);
                    }
                }
            }
        }

        // 获取去重后所有点在level0的入度
        std::cout << "call indegree level0" << std::endl;
        std::vector<linklistsizeint> in_degrees(max_elements_);
        // for (int i = 0; i < merge_graph.size(); i++) {
        //     auto& internal_neighbours = merge_graph[i].internal_neighbours_;
        //     for (auto& internal_neighbours_id: internal_neighbours[0]) {
        //         in_degrees[internal_neighbours_id]++;
        //     }
        // }

        std::cout << "fill edges" << std::endl;
        // 填充邻居
        // std::random_device rng;
        // std::mt19937 urng(rng());
        // size_links_per_element_ = maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);  // 高层的 level 也变成 2*M 实时看
        size_links_per_element_ = maxM_ * sizeof(tableint) + sizeof(linklistsizeint);

        std::cout << "enter point node: " << enterpoint_node_ << std::endl;
        size_t enterpoint_max_level = 0;
        for (int i = 0; i < merge_graph.size(); i++) {
            auto internal_label = merge_graph[i].internal_label_;
            auto cur_max_level = merge_graph[i].max_level_;
            element_levels_[i] = cur_max_level;
            if (cur_max_level > enterpoint_max_level) {
                enterpoint_max_level = cur_max_level;
                enterpoint_node_ = internal_label;
            }

            if (cur_max_level > 0) {
                linkLists_[internal_label] = (char*)malloc(size_links_per_element_*cur_max_level+1);
                memset(linkLists_[internal_label], 0, size_links_per_element_*cur_max_level+1);
            }

            auto& internal_neighbours = merge_graph[i].internal_neighbours_;
            for (int level = 0; level <= cur_max_level; level++) {
                auto& internal_neighbours_level = internal_neighbours[level];
                mergeSelectNeighbors(i, internal_neighbours_level, level, in_degrees);
                // std::shuffle(internal_neighbours_level.begin(), internal_neighbours_level.end(), urng);

                if (level == 0) {
                    size_t num_neighbours = internal_neighbours_level.size() > maxM0_ ? maxM0_ : internal_neighbours_level.size();

                    linklistsizeint* ll_cur = get_linklist0(internal_label);
                    setListCount(ll_cur, num_neighbours);
                    tableint* data = (tableint*)(ll_cur+1);
                    for (int i = 0; i < num_neighbours; i++) {
                        data[i] = internal_neighbours_level[i];
                    }
                } else {
                    size_t num_neighbours = internal_neighbours_level.size() > maxM_ ? maxM_ : internal_neighbours_level.size();

                    linklistsizeint* ll_cur = get_linklist(internal_label, level);
                    setListCount(ll_cur, num_neighbours);
                    tableint* data = (tableint*)(ll_cur+1);
                    for (int i = 0; i < num_neighbours; i++) {
                        data[i] = internal_neighbours_level[i];
                    }
                }
            }
        }
        std::cout << "enter point node: " << enterpoint_node_ << std::endl;

        // // 增加反向连边，
        // // 方案 1. 在构图结束后，从 id=0 开始进行反向连接
        // // 方案 2. 在构图结束后，根据入度进行排列，入度小的首先进行选择
        // // 方案 3. 在构图前根据入度进行排序，然后选择 
        // std::cout << "add reverse edges" << std::endl;
        // for (int i = 0; i < merge_graph.size(); i++) {
        //     auto max_level = element_levels_[i];
        //     for (int level = 0; level <= max_level; level++) {
        //         linklistsizeint* ll_cur = get_linklist_by_level(i, level);
        //         tableint* data_neighbours = (tableint*)(ll_cur+1);
        //         for (linklistsizeint j = 0; j < *ll_cur; j++) {
        //             tableint internal_id_neighbour = data_neighbours[j];
        //             linklistsizeint* ll_cur_nei = get_linklist_by_level(internal_id_neighbour, level);
        //             tableint* data_nei = (tableint*)(ll_cur_nei+1);
        //
        //             auto m = level == 0 ? maxM0_ : maxM_;
        //             if (*ll_cur_nei < m) {
        //                 bool found = false;
        //                 for (linklistsizeint k = 0; k < *ll_cur_nei; k++) {
        //                     if (data_nei[k] == i) {
        //                         found = true;
        //                         break;
        //                     }
        //                 }
        //
        //                 if (!found) {
        //                     data_nei[*ll_cur_nei] = i;
        //                     *ll_cur_nei = *ll_cur_nei + 1;
        //                 }
        //             }
        //         }
        //     }
        // }

        // std::cout << "cal in degree" << std::endl;
        // std::vector<int> indices(max_elements_);
        // std::vector<int> in_degree(max_elements_, 0);
        // // 计算level0入度，并排序，从入度最小的开始选择
        // for (tableint i = 0; i < max_elements_; i++) {
        //     indices[i] = i;
        //     linklistsizeint* ll_cur = get_linklist0(i);
        //
        //     tableint num_neighbours = *ll_cur;
        //     tableint* data_neighbours = (tableint*)(ll_cur+1);
        //     for (tableint j = 0; j < num_neighbours; j++) {
        //         tableint internal_id_neighbour = data_neighbours[j];
        //         in_degree[internal_id_neighbour]++;
        //     }
        // }
        //
        // std::cout << "sort in degree" << std::endl;
        // std::sort(indices.begin(), indices.end(),
        //       [&in_degree](int i1, int i2) { return in_degree[i1] < in_degree[i2]; });
        //
        // std::cout << "reconnect" << std::endl;
        // for (tableint i = 0; i < max_elements_; i++) {
        //     for (int level = 0; level <= element_levels_[indices[i]]; level++) {
        //         std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
        //         linklistsizeint* ll_cur = get_linklist_by_level(indices[i], level);
        //
        //         tableint num_neighbours = *ll_cur;
        //         tableint* data_neighbours = (tableint*)(ll_cur+1);
        //         for (tableint j = 0; j < num_neighbours; j++) {
        //             tableint internal_id_neighbour = data_neighbours[j];
        //             top_candidates.push(std::make_pair(0.0, internal_id_neighbour));
        //         }
        //
        //         mutuallyConnectNewElement(nullptr, indices[i], top_candidates, level, false, false);
        //     }
        // }
    }

    linklistsizeint *get_linklist_by_level(tableint internal_id, int level) const {
        if (level == 0) {
            return get_linklist0(internal_id);
        } else {
            return get_linklist(internal_id, level);
        }
    }


    void mergeSelectNeighbors(tableint home, std::vector<tableint>& internal_neighbours, int level, std::vector<linklistsizeint>& in_degree) {


        // 方案 1：去重后随机选择
        // 提前去重
        {
            int current_m = level == 0 ? maxM0_ : maxM_;
            std::sort(internal_neighbours.begin(), internal_neighbours.end());
            auto it = std::unique(internal_neighbours.begin(), internal_neighbours.end());
            internal_neighbours.erase(it, internal_neighbours.end());
            if (internal_neighbours.size() <= current_m) {
                return;
            }
            std::random_device rng;
            std::mt19937 urng(rng());
            std::shuffle(internal_neighbours.begin(), internal_neighbours.end(), urng);
            internal_neighbours.resize(current_m);
        }

        // // 方案 2：按照 id 排序，重复多的放到最前面
        // int current_m = level == 0 ? maxM0_ : maxM_;
        // if ((level == 0 && internal_neighbours.size() < maxM0_) ||
        //     (level > 0 && internal_neighbours.size() < maxM_)) {
        //     // 数量少直接去重
        //     std::sort(internal_neighbours.begin(), internal_neighbours.end());
        //     auto it = std::unique(internal_neighbours.begin(), internal_neighbours.end());
        //     internal_neighbours.erase(it, internal_neighbours.end());
        //     return;
        // }
        //
        // // 统计重复数量
        // std::unordered_map<tableint, int> count_map;
        // for (auto& internal_id : internal_neighbours) {
        //     count_map[internal_id]++;
        // }
        // // 根据数量排序
        // std::vector<std::pair<tableint, int>> count_neighbours(count_map.begin(), count_map.end());
        // std::sort(count_neighbours.begin(), count_neighbours.end(), [](const auto &left, const auto &right) {
        //     return left.second > right.second;
        // });
        //
        // std::vector<tableint>().swap(internal_neighbours);  // 清空并释放内存
        // internal_neighbours.push_back(count_neighbours[0].first);  // 重复最多的放到最前面
        // for (int i = 1; i < count_neighbours.size() && internal_neighbours.size() < current_m; i++) {
        //     if (count_neighbours[i].first != count_neighbours[i-1].first) {
        //         internal_neighbours.push_back(count_neighbours[i].first);
        //     }
        // }

        // // 方案 3，按照邻居的入度来删减, 邻居入度少的优先留下
        // // 需要整个图提前去重
        // {
        //     int current_m = level == 0 ? maxM0_ : maxM_;
        //     if ((level == 0 && internal_neighbours.size() <= maxM0_) ||
        //         (level > 0 && internal_neighbours.size() <= maxM_)) {
        //         // 数量少直接返回
        //         return;
        //     }
        //
        //     // 上层的采用随机选择的方案
        //     if (level > 0) {
        //         std::random_device rng;
        //         std::mt19937 urng(rng());
        //         std::shuffle(internal_neighbours.begin(), internal_neighbours.end(), urng);
        //         internal_neighbours.resize(maxM_);
        //         return ;
        //     }
        //
        //     // level0 根据入度排序进行选择
        //     std::vector<std::pair<tableint, int>> neighbour_with_indegree(internal_neighbours.size());
        //     for (int i = 0; i < internal_neighbours.size(); i++) {
        //         neighbour_with_indegree[i] = std::make_pair(internal_neighbours[i], in_degree[internal_neighbours[i]]);
        //     }
        //     std::sort(neighbour_with_indegree.begin(), neighbour_with_indegree.end(), [](const auto &left, const auto &right) {
        //             return left.second > right.second;
        //     });
        //     internal_neighbours.resize(maxM0_);
        //     // std::copy(neighbour_with_indegree.begin(), neighbour_with_indegree.begin() + maxM0_, internal_neighbours.begin());
        //     // std::vector<tableint>().swap(internal_neighbours);  // 清空并释放内存
        //     for (int i = 0; i < neighbour_with_indegree.size(); i++) {
        //         if (i < maxM0_) {
        //             internal_neighbours[i] = neighbour_with_indegree[i].first;
        //         } else {
        //             in_degree[neighbour_with_indegree[i].first]--;
        //         }
        //     }
        // }

//         // 方案 3：启发式选择
//         if (internal_neighbours.size() >= 32 ) {
//             return;
//         }
//         std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
// // #pragma omp parallel for
//         for (int i = 0; i < internal_neighbours.size(); i++) {
//             tableint neighbour_internal_id = internal_neighbours[i];
//             dist_t distance = fstdistfunc_(
//                     getDataByInternalId(home), 
//                     getDataByInternalId(neighbour_internal_id), dist_func_param_, scale2_);
//             top_candidates.push(std::make_pair(distance, neighbour_internal_id));
//         }
//
//         getNeighborsByHeuristic2(top_candidates, level == 0 ? maxM0_ : maxM_);
//         std::vector<tableint>().swap(internal_neighbours);
//         while (!top_candidates.empty()) {
//             internal_neighbours.push_back(top_candidates.top().second);
//             top_candidates.pop();
//         }

        // 方案 4
        // 增加反向连边
    }

    // // 计算合并前所有层点的入度
    // void calInDegrees0(std::vector<linklistsizeint>& in_degree) {
    //     for (int i = 0; i < max_elements_; i++) {
    //         linklistsizeint* ll_cur = get_linklist0(i);
    //
    //         linklistsizeint num_neighbours = *ll_cur;
    //         tableint* data_neighbours = (tableint*)(ll_cur+1);
    //         for (tableint j = 0; j < num_neighbours; j++) {
    //             tableint internal_id_neighbour = data_neighbours[j];
    //             in_degree[internal_id_neighbour]++;
    //         }
    //     }
    // }

    void loadCodeBooks(const std::vector<std::vector<float>>& code_books) {
        pq_centroids_ = code_books;
        hnswlib::codebooks = code_books;
    }

    void loadResiduals(const std::vector<float>& residuals) {
        pq_residuals_ = residuals;
    }

    // void loadQueryResiduals(const std::vector<float>& query_residuals) {
    //     query_residuals_ = query_residuals;
    // }

    void loadPqIndex(const std::vector<std::vector<uint8_t>>& pq_codes) {
        for (int i = 0; i < max_elements_; i++) {
            char* data = getDataByInternalId(i);
            memcpy(data, pq_codes[i].data(), pq_codes[i].size());
        }
    }

    float L2(const float* d1, const float *d2, int dim) {
        float r = 0.0;
        for (int i = 0; i < dim; i++) {
            float r1 = *(d1+i)-*(d2+i);
            r += r1*r1;
        }

        return r;
    }

    void calDistLookUpTable(int M, int Ks, int dsub) {
      hnswlib::dist_lookup.resize(M);
      for (int i = 0; i < hnswlib::dist_lookup.size(); i++) {
          hnswlib::dist_lookup[i].resize(Ks * (Ks + 1) / 2);
      }

      for (int i = 0; i < M; i++) {
        const float *data = pq_centroids_[i].data();

        for (int j = 0; j < Ks; j++) {
          for (int k = 0; k < Ks; k++) {
            if (j < k)
              break;

            float d = L2(data + j * dsub, data + k * dsub, dsub);
            size_t idx = j * (j + 1) / 2 + k;
            dist_lookup[i][idx] = d;
          }
        }
      }
    }

    float calMax() {
        float max_ = 0.0;
        float per = 0.9;
        int64_t top_n = (1.0-per)*max_elements_;

        std::priority_queue<float, std::vector<float>, std::greater<float>> top_candidates;
        size_t dim = *((size_t *) dist_func_param_);
        for (size_t i = 0; i < max_elements_; i++) {
            char* data_ptrv = getDataByInternalId(i);
            float* data_ptr = (float*) data_ptrv;
            for (size_t i = 0; i < dim; i++) {
                float abs_val = data_ptr[i] > 0 ? data_ptr[i] : -data_ptr[i];
                if (top_candidates.size() < top_n) {
                    top_candidates.push(abs_val);
                } 
                if (top_candidates.size() == top_n) {
                    float min_val = top_candidates.top();
                    if (abs_val > min_val) {
                        top_candidates.pop();
                        top_candidates.push(abs_val);
                    }
                }
            }
        }

        float max_val = top_candidates.top();

        while (max_val <= 0 && !top_candidates.empty()) {
            max_val = top_candidates.top();
            top_candidates.pop();
        }

        return max_val;
    }

    float getScale() {
        return scale_;
    }

    void convertVectorInplace(char* data_ptr) {
        size_t dim = *((size_t *) dist_func_param_);
        std::vector<int8_t> tmp_v(dim);

        float* data_ptrf = (float*) data_ptr;
        for (size_t i = 0; i < dim; i++) {
            float scaled_x = data_ptrf[i] * scale_;
            if (scaled_x > 127) {
                scaled_x = 127;
            } else if (scaled_x <= -128) {
                scaled_x = -128;
            }

            tmp_v[i] = (int8_t)(scaled_x);
        }

        memcpy(data_ptr, tmp_v.data(), dim);
    }

    void sq8() {
        float max_val = calMax();
        scale_ = 127 / max_val;
        scale2_ = scale_ * scale_;
#pragma omp parallel for
        for (int i = 0; i < max_elements_; i++) {
            char* data_ptrv = getDataByInternalId(i);
            convertVectorInplace(data_ptrv );
        }
    }

    template<typename data_t>
    std::vector<data_t> getDataByLabel(labeltype label) const {
        // lock all operations with element by label
        std::unique_lock <std::mutex> lock_label(getLabelOpMutex(label));
        
        std::unique_lock <std::mutex> lock_table(label_lookup_lock);
        auto search = label_lookup_.find(label);
        if (search == label_lookup_.end() || isMarkedDeleted(search->second)) {
            throw std::runtime_error("Label not found");
        }
        tableint internalId = search->second;
        lock_table.unlock();

        char* data_ptrv = getDataByInternalId(internalId);
        size_t dim = *((size_t *) dist_func_param_);
        std::vector<data_t> data;
        data_t* data_ptr = (data_t*) data_ptrv;
        for (size_t i = 0; i < dim; i++) {
            data.push_back(*data_ptr);
            data_ptr += 1;
        }
        return data;
    }

    std::vector<float> getDataByLabelFloat(labeltype label) const {
        // lock all operations with element by label
        std::unique_lock <std::mutex> lock_label(getLabelOpMutex(label));
        
        std::unique_lock <std::mutex> lock_table(label_lookup_lock);
        auto search = label_lookup_.find(label);
        if (search == label_lookup_.end() || isMarkedDeleted(search->second)) {
            throw std::runtime_error("Label not found");
        }
        tableint internalId = search->second;
        lock_table.unlock();

        char* data_ptrv = getDataByInternalId(internalId);
        size_t dim = *((size_t *) dist_func_param_);
        std::vector<float> data;
        float* data_ptr = (float*) data_ptrv;
        for (size_t i = 0; i < dim; i++) {
            data.push_back(*data_ptr);
            data_ptr += 1;
        }
        return data;
    }


    /*
    * Marks an element with the given label deleted, does NOT really change the current graph.
    */
    void markDelete(labeltype label) {
        // lock all operations with element by label
        std::unique_lock <std::mutex> lock_label(getLabelOpMutex(label));

        std::unique_lock <std::mutex> lock_table(label_lookup_lock);
        auto search = label_lookup_.find(label);
        if (search == label_lookup_.end()) {
            throw std::runtime_error("Label not found");
        }
        tableint internalId = search->second;
        lock_table.unlock();

        markDeletedInternal(internalId);
    }


    /*
    * Uses the last 16 bits of the memory for the linked list size to store the mark,
    * whereas maxM0_ has to be limited to the lower 16 bits, however, still large enough in almost all cases.
    */
    void markDeletedInternal(tableint internalId) {
        assert(internalId < cur_element_count);
        if (!isMarkedDeleted(internalId)) {
            unsigned char *ll_cur = ((unsigned char *)get_linklist0(internalId))+2;
            *ll_cur |= DELETE_MARK;
            num_deleted_ += 1;
            if (allow_replace_deleted_) {
                std::unique_lock <std::mutex> lock_deleted_elements(deleted_elements_lock);
                deleted_elements.insert(internalId);
            }
        } else {
            throw std::runtime_error("The requested to delete element is already deleted");
        }
    }


    /*
    * Removes the deleted mark of the node, does NOT really change the current graph.
    * 
    * Note: the method is not safe to use when replacement of deleted elements is enabled,
    *  because elements marked as deleted can be completely removed by addPoint
    */
    void unmarkDelete(labeltype label) {
        // lock all operations with element by label
        std::unique_lock <std::mutex> lock_label(getLabelOpMutex(label));

        std::unique_lock <std::mutex> lock_table(label_lookup_lock);
        auto search = label_lookup_.find(label);
        if (search == label_lookup_.end()) {
            throw std::runtime_error("Label not found");
        }
        tableint internalId = search->second;
        lock_table.unlock();

        unmarkDeletedInternal(internalId);
    }



    /*
    * Remove the deleted mark of the node.
    */
    void unmarkDeletedInternal(tableint internalId) {
        assert(internalId < cur_element_count);
        if (isMarkedDeleted(internalId)) {
            unsigned char *ll_cur = ((unsigned char *)get_linklist0(internalId)) + 2;
            *ll_cur &= ~DELETE_MARK;
            num_deleted_ -= 1;
            if (allow_replace_deleted_) {
                std::unique_lock <std::mutex> lock_deleted_elements(deleted_elements_lock);
                deleted_elements.erase(internalId);
            }
        } else {
            throw std::runtime_error("The requested to undelete element is not deleted");
        }
    }


    /*
    * Checks the first 16 bits of the memory to see if the element is marked deleted.
    */
    bool isMarkedDeleted(tableint internalId) const {
        unsigned char *ll_cur = ((unsigned char*)get_linklist0(internalId)) + 2;
        return *ll_cur & DELETE_MARK;
    }


    unsigned short int getListCount(linklistsizeint * ptr) const {
        return *((unsigned short int *)ptr);
    }


    void setListCount(linklistsizeint * ptr, unsigned short int size) const {
        *((unsigned short int*)(ptr))=*((unsigned short int *)&size);
    }


    /*
    * Adds point. Updates the point if it is already in the index.
    * If replacement of deleted elements is enabled: replaces previously deleted point if any, updating it with new point
    */
    void addPoint(const void *data_point, labeltype label, bool replace_deleted = false) {
        if ((allow_replace_deleted_ == false) && (replace_deleted == true)) {
            throw std::runtime_error("Replacement of deleted elements is disabled in constructor");
        }

        // lock all operations with element by label
        std::unique_lock <std::mutex> lock_label(getLabelOpMutex(label));
        if (!replace_deleted) {
            addPoint(data_point, label, -1);
            return;
        }
        // check if there is vacant place
        tableint internal_id_replaced;
        std::unique_lock <std::mutex> lock_deleted_elements(deleted_elements_lock);
        bool is_vacant_place = !deleted_elements.empty();
        if (is_vacant_place) {
            internal_id_replaced = *deleted_elements.begin();
            deleted_elements.erase(internal_id_replaced);
        }
        lock_deleted_elements.unlock();

        // if there is no vacant place then add or update point
        // else add point to vacant place
        if (!is_vacant_place) {
            addPoint(data_point, label, -1);
        } else {
            // we assume that there are no concurrent operations on deleted element
            labeltype label_replaced = getExternalLabel(internal_id_replaced);
            setExternalLabel(internal_id_replaced, label);

            std::unique_lock <std::mutex> lock_table(label_lookup_lock);
            label_lookup_.erase(label_replaced);
            label_lookup_[label] = internal_id_replaced;
            lock_table.unlock();

            unmarkDeletedInternal(internal_id_replaced);
            updatePoint(data_point, internal_id_replaced, 1.0);
        }
    }


    void updatePoint(const void *dataPoint, tableint internalId, float updateNeighborProbability) {
        // update the feature vector associated with existing point with new vector
        memcpy(getDataByInternalId(internalId), dataPoint, data_size_);

        int maxLevelCopy = maxlevel_;
        tableint entryPointCopy = enterpoint_node_;
        // If point to be updated is entry point and graph just contains single element then just return.
        if (entryPointCopy == internalId && cur_element_count == 1)
            return;

        int elemLevel = element_levels_[internalId];
        std::uniform_real_distribution<float> distribution(0.0, 1.0);
        for (int layer = 0; layer <= elemLevel; layer++) {
            std::unordered_set<tableint> sCand;
            std::unordered_set<tableint> sNeigh;
            std::vector<tableint> listOneHop = getConnectionsWithLock(internalId, layer);
            if (listOneHop.size() == 0)
                continue;

            sCand.insert(internalId);

            for (auto&& elOneHop : listOneHop) {
                sCand.insert(elOneHop);

                if (distribution(update_probability_generator_) > updateNeighborProbability)
                    continue;

                sNeigh.insert(elOneHop);

                std::vector<tableint> listTwoHop = getConnectionsWithLock(elOneHop, layer);
                for (auto&& elTwoHop : listTwoHop) {
                    sCand.insert(elTwoHop);
                }
            }

            for (auto&& neigh : sNeigh) {
                // if (neigh == internalId)
                //     continue;

                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidates;
                size_t size = sCand.find(neigh) == sCand.end() ? sCand.size() : sCand.size() - 1;  // sCand guaranteed to have size >= 1
                size_t elementsToKeep = std::min(ef_construction_, size);
                for (auto&& cand : sCand) {
                    if (cand == neigh)
                        continue;

                    dist_t distance = fstdistfunc_(getDataByInternalId(neigh), getDataByInternalId(cand), dist_func_param_, scale2_);
                    if (candidates.size() < elementsToKeep) {
                        candidates.emplace(distance, cand);
                    } else {
                        if (distance < candidates.top().first) {
                            candidates.pop();
                            candidates.emplace(distance, cand);
                        }
                    }
                }

                // Retrieve neighbours using heuristic and set connections.
                getNeighborsByHeuristic2(candidates, layer == 0 ? maxM0_ : maxM_);

                {
                    std::unique_lock <std::mutex> lock(link_list_locks_[neigh]);
                    linklistsizeint *ll_cur;
                    ll_cur = get_linklist_at_level(neigh, layer);
                    size_t candSize = candidates.size();
                    setListCount(ll_cur, candSize);
                    tableint *data = (tableint *) (ll_cur + 1);
                    for (size_t idx = 0; idx < candSize; idx++) {
                        data[idx] = candidates.top().second;
                        candidates.pop();
                    }
                }
            }
        }

        repairConnectionsForUpdate(dataPoint, entryPointCopy, internalId, elemLevel, maxLevelCopy);
    }


    void repairConnectionsForUpdate(
        const void *dataPoint,
        tableint entryPointInternalId,
        tableint dataPointInternalId,
        int dataPointLevel,
        int maxLevel) {
        tableint currObj = entryPointInternalId;
        if (dataPointLevel < maxLevel) {
            dist_t curdist = fstdistfunc_(dataPoint, getDataByInternalId(currObj), dist_func_param_, scale2_);
            for (int level = maxLevel; level > dataPointLevel; level--) {
                bool changed = true;
                while (changed) {
                    changed = false;
                    unsigned int *data;
                    std::unique_lock <std::mutex> lock(link_list_locks_[currObj]);
                    data = get_linklist_at_level(currObj, level);
                    int size = getListCount(data);
                    tableint *datal = (tableint *) (data + 1);
#ifdef USE_SSE
                    _mm_prefetch(getDataByInternalId(*datal), _MM_HINT_T0);
#endif
                    for (int i = 0; i < size; i++) {
#ifdef USE_SSE
                        _mm_prefetch(getDataByInternalId(*(datal + i + 1)), _MM_HINT_T0);
#endif
                        tableint cand = datal[i];
                        dist_t d = fstdistfunc_(dataPoint, getDataByInternalId(cand), dist_func_param_, scale2_);
                        if (d < curdist) {
                            curdist = d;
                            currObj = cand;
                            changed = true;
                        }
                    }
                }
            }
        }

        if (dataPointLevel > maxLevel)
            throw std::runtime_error("Level of item to be updated cannot be bigger than max level");

        for (int level = dataPointLevel; level >= 0; level--) {
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> topCandidates = searchBaseLayer(
                    currObj, dataPoint, level);

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> filteredTopCandidates;
            while (topCandidates.size() > 0) {
                if (topCandidates.top().second != dataPointInternalId)
                    filteredTopCandidates.push(topCandidates.top());

                topCandidates.pop();
            }

            // Since element_levels_ is being used to get `dataPointLevel`, there could be cases where `topCandidates` could just contains entry point itself.
            // To prevent self loops, the `topCandidates` is filtered and thus can be empty.
            if (filteredTopCandidates.size() > 0) {
                bool epDeleted = isMarkedDeleted(entryPointInternalId);
                if (epDeleted) {
                    filteredTopCandidates.emplace(fstdistfunc_(dataPoint, getDataByInternalId(entryPointInternalId), dist_func_param_, scale2_), entryPointInternalId);
                    if (filteredTopCandidates.size() > ef_construction_)
                        filteredTopCandidates.pop();
                }

                currObj = mutuallyConnectNewElement(dataPoint, dataPointInternalId, filteredTopCandidates, level, true);
            }
        }
    }


    std::vector<tableint> getConnectionsWithLock(tableint internalId, int level) {
        std::unique_lock <std::mutex> lock(link_list_locks_[internalId]);
        unsigned int *data = get_linklist_at_level(internalId, level);
        int size = getListCount(data);
        std::vector<tableint> result(size);
        tableint *ll = (tableint *) (data + 1);
        memcpy(result.data(), ll, size * sizeof(tableint));
        return result;
    }


    tableint addPoint(const void *data_point, labeltype label, int level) {
        tableint cur_c = 0;
        {
            // Checking if the element with the same label already exists
            // if so, updating it *instead* of creating a new element.
            std::unique_lock <std::mutex> lock_table(label_lookup_lock);
            auto search = label_lookup_.find(label);
            if (search != label_lookup_.end()) {
                tableint existingInternalId = search->second;
                if (allow_replace_deleted_) {
                    if (isMarkedDeleted(existingInternalId)) {
                        throw std::runtime_error("Can't use addPoint to update deleted elements if replacement of deleted elements is enabled.");
                    }
                }
                lock_table.unlock();

                if (isMarkedDeleted(existingInternalId)) {
                    unmarkDeletedInternal(existingInternalId);
                }
                updatePoint(data_point, existingInternalId, 1.0);

                return existingInternalId;
            }

            if (cur_element_count >= max_elements_) {
                throw std::runtime_error("The number of elements exceeds the specified limit");
            }

            cur_c = cur_element_count;
            cur_element_count++;
            label_lookup_[label] = cur_c;
        }

        std::unique_lock <std::mutex> lock_el(link_list_locks_[cur_c]);
        int curlevel = getRandomLevel(mult_);
        // int curlevel = getRandomLevel(revSize_);
        if (level > 0)
            curlevel = level;

        element_levels_[cur_c] = curlevel;  // 设置当前点所属的 level，最高 level

        std::unique_lock <std::mutex> templock(global);
        int maxlevelcopy = maxlevel_;
        if (curlevel <= maxlevelcopy)
            templock.unlock();
        tableint currObj = enterpoint_node_;
        tableint enterpoint_copy = enterpoint_node_;  // 获取进入点

        memset(data_level0_memory_ + cur_c * size_data_per_element_ + offsetLevel0_, 0, size_data_per_element_);

        // Initialisation of the data and label
        memcpy(getExternalLabeLp(cur_c), &label, sizeof(labeltype)); // level0 写入外部 id
        memcpy(getDataByInternalId(cur_c), data_point, data_size_);  // level0 写入数据

        if (curlevel) {  // 如果当前点不是在第 0 层
            linkLists_[cur_c] = (char *) malloc(size_links_per_element_ * curlevel + 1); // 为这个点分配curlevel个层，每个层都有 M 个邻居
            if (linkLists_[cur_c] == nullptr)
                throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
            memset(linkLists_[cur_c], 0, size_links_per_element_ * curlevel + 1);
        }

        if ((signed)currObj != -1) {
            if (curlevel < maxlevelcopy) { // 寻找到当前层的进入点
                dist_t curdist = fstdistfunc_(data_point, getDataByInternalId(currObj), dist_func_param_, scale2_);
                for (int level = maxlevelcopy; level > curlevel; level--) {
                    bool changed = true;
                    while (changed) {
                        changed = false;
                        unsigned int *data;
                        std::unique_lock <std::mutex> lock(link_list_locks_[currObj]);
                        data = get_linklist(currObj, level);
                        int size = getListCount(data);

                        tableint *datal = (tableint *) (data + 1);
                        for (int i = 0; i < size; i++) {
                            tableint cand = datal[i];
                            if (cand < 0 || cand > max_elements_)
                                throw std::runtime_error("cand error");
                            dist_t d = fstdistfunc_(data_point, getDataByInternalId(cand), dist_func_param_, scale2_);
                            if (d < curdist) {
                                curdist = d;
                                currObj = cand;
                                changed = true;
                            }
                        }
                    }
                }
            }

            bool epDeleted = isMarkedDeleted(enterpoint_copy);
            for (int level = std::min(curlevel, maxlevelcopy); level >= 0; level--) {  // 从当前层开始进入，为每一层构图，知道 level0
                if (level > maxlevelcopy || level < 0)  // possible?
                    throw std::runtime_error("Level error");

                // 从当前层的进入点开始寻找最近的邻居
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates = searchBaseLayer(
                        currObj, data_point, level);
                if (epDeleted) {
                    top_candidates.emplace(fstdistfunc_(data_point, getDataByInternalId(enterpoint_copy), dist_func_param_, scale2_), enterpoint_copy);
                    if (top_candidates.size() > ef_construction_)
                        top_candidates.pop();
                }
                // 在当前层构图，同时获得下一层的进入点 
                currObj = mutuallyConnectNewElement(data_point, cur_c, top_candidates, level, false); // 为 level 层创建连接
            }
        } else {
            // Do nothing for the first element
            enterpoint_node_ = 0;
            maxlevel_ = curlevel;
        }

        // Releasing lock for the maximum level
        if (curlevel > maxlevelcopy) {
            enterpoint_node_ = cur_c;
            maxlevel_ = curlevel;
        }
        return cur_c;
    }


    std::priority_queue<std::pair<dist_t, labeltype >>
    searchKnn(const void *query_data, size_t k, float q_residual, BaseFilterFunctor* isIdAllowed = nullptr) const {
        std::priority_queue<std::pair<dist_t, labeltype >> result;
        if (cur_element_count == 0) return result;

        tableint currObj = enterpoint_node_;
        dist_t curdist = fstdistfunc_(query_data, getDataByInternalId(enterpoint_node_), dist_func_param_, scale2_);
        // add residuals
        // curdist += q_residual;
        // curdist += pq_residuals_[*getExternalLabeLp(enterpoint_node_)];

        for (int level = maxlevel_; level > 0; level--) {
            bool changed = true;
            while (changed) {
                changed = false;
                unsigned int *data;

                data = (unsigned int *) get_linklist(currObj, level);
                int size = getListCount(data);
                metric_hops++;
                metric_distance_computations+=size;

                tableint *datal = (tableint *) (data + 1);
                for (int i = 0; i < size; i++) {
                    tableint cand = datal[i];
                    if (cand < 0 || cand > max_elements_)
                        throw std::runtime_error("cand error");
                    dist_t d = fstdistfunc_(query_data, getDataByInternalId(cand), dist_func_param_, scale2_);
                    // add residuals
                    // d += q_residual;
                    // d += pq_residuals_[*getExternalLabeLp(cand)];

                    if (d < curdist) {
                        curdist = d;
                        currObj = cand;
                        changed = true;
                    }
                }
            }
        }

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
        bool bare_bone_search = !num_deleted_ && !isIdAllowed;
        if (bare_bone_search) {
            top_candidates = searchBaseLayerST<true, true>( // collect_metrics
                    currObj, query_data, std::max(ef_, k), q_residual, isIdAllowed);
        } else {
            top_candidates = searchBaseLayerST<false>(
                    currObj, query_data, std::max(ef_, k), q_residual, isIdAllowed);
        }

        while (top_candidates.size() > k) {
            top_candidates.pop();
        }
        while (top_candidates.size() > 0) {
            std::pair<dist_t, tableint> rez = top_candidates.top();
            result.push(std::pair<dist_t, labeltype>(rez.first, getExternalLabel(rez.second)));
            top_candidates.pop();
        }
        return result;
    }


    std::vector<std::pair<dist_t, labeltype >>
    searchStopConditionClosest(
        const void *query_data,
        BaseSearchStopCondition<dist_t>& stop_condition,
        BaseFilterFunctor* isIdAllowed = nullptr) const {
        std::vector<std::pair<dist_t, labeltype >> result;
        if (cur_element_count == 0) return result;

        tableint currObj = enterpoint_node_;
        dist_t curdist = fstdistfunc_(query_data, getDataByInternalId(enterpoint_node_), dist_func_param_, scale2_);

        for (int level = maxlevel_; level > 0; level--) {
            bool changed = true;
            while (changed) {
                changed = false;
                unsigned int *data;

                data = (unsigned int *) get_linklist(currObj, level);
                int size = getListCount(data);
                metric_hops++;
                metric_distance_computations+=size;

                tableint *datal = (tableint *) (data + 1);
                for (int i = 0; i < size; i++) {
                    tableint cand = datal[i];
                    if (cand < 0 || cand > max_elements_)
                        throw std::runtime_error("cand error");
                    dist_t d = fstdistfunc_(query_data, getDataByInternalId(cand), dist_func_param_, scale2_);

                    if (d < curdist) {
                        curdist = d;
                        currObj = cand;
                        changed = true;
                    }
                }
            }
        }

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
        top_candidates = searchBaseLayerST<false>(currObj, query_data, 0, isIdAllowed, &stop_condition);

        size_t sz = top_candidates.size();
        result.resize(sz);
        while (!top_candidates.empty()) {
            result[--sz] = top_candidates.top();
            top_candidates.pop();
        }

        stop_condition.filter_results(result);

        return result;
    }


    void checkIntegrity() {
        int connections_checked = 0;
        std::vector <int > inbound_connections_num(cur_element_count, 0);
        for (int i = 0; i < cur_element_count; i++) {
            for (int l = 0; l <= element_levels_[i]; l++) {
                linklistsizeint *ll_cur = get_linklist_at_level(i, l);
                int size = getListCount(ll_cur);
                tableint *data = (tableint *) (ll_cur + 1);
                std::unordered_set<tableint> s;
                for (int j = 0; j < size; j++) {
                    assert(data[j] < cur_element_count);
                    assert(data[j] != i);
                    inbound_connections_num[data[j]]++;
                    s.insert(data[j]);
                    connections_checked++;
                }
                assert(s.size() == size);
            }
        }
        if (cur_element_count > 1) {
            int min1 = inbound_connections_num[0], max1 = inbound_connections_num[0];
            for (int i=0; i < cur_element_count; i++) {
                assert(inbound_connections_num[i] > 0);
                min1 = std::min(inbound_connections_num[i], min1);
                max1 = std::max(inbound_connections_num[i], max1);
            }
            std::cout << "Min inbound: " << min1 << ", Max inbound:" << max1 << "\n";
        }
        std::cout << "integrity ok, checked " << connections_checked << " connections\n";
    }
};
}  // namespace hnswlib
   //
