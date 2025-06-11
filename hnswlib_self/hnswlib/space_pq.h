#pragma once
#include "./hnswlib.h"

namespace hnswlib {

std::vector<std::vector<float>> codebooks;
std::vector<std::vector<float>> dist_lookup;

static float
sdc_pq_distance(const void *pVect1v, const void *pVect2v, const void *qty_ptr, float t) {
  uint8_t *pv1 = (uint8_t *)pVect1v;
  uint8_t *pv2 = (uint8_t *)pVect2v;

  size_t qty = *((size_t *)qty_ptr);
  float res = 0;

  for (size_t i = 0; i < qty; ++i) {
    uint8_t idx1 = pv1[i];
    uint8_t idx2 = pv2[i];
    if (idx1 < idx2) {
        size_t idx = idx2 * (idx2 + 1) / 2 + idx1;
        res += dist_lookup[i][idx];
    } else {
        size_t idx = idx1 * (idx1 + 1) / 2 + idx2;
        res += dist_lookup[i][idx];
    }
  }
  return res;
}

static float adc_pq_distance(const void *pVect1v, const void *pVect2v,
                             const void *qty_ptr, float t) {
  float *pVect1 = (float *)pVect1v;
  uint8_t *pVect2 = (uint8_t *)pVect2v;
  size_t M = *((size_t *)qty_ptr);

  int dim = 128;
  int dsub = dim / M;

  float res = 0;
  for (size_t i = 0; i < M; i++) {
    uint8_t code_books_idx = pVect2[i];
    float* centroids = codebooks[i].data()+code_books_idx*dsub;

    for (size_t j = 0; j < dsub; j++) {
        float d1 = pVect1[i * dsub + j];
        float d2 = centroids[j];
        res += (d1 - d2) * (d1 - d2);
    }
  }

  return res;
}

class PqSpace : public SpaceInterface<float> {
    DISTFUNC<float> fstdistfunc_;
    size_t data_size_;
    size_t dim_;

 public:
    PqSpace(size_t dim) {
        fstdistfunc_ = adc_pq_distance;
        dim_ = dim;
        data_size_ = dim * sizeof(uint8_t);
    }

    size_t get_data_size() {
        return data_size_;
    }

    DISTFUNC<float> get_dist_func() {
        return fstdistfunc_;
    }

    void *get_dist_func_param() {
        return &dim_;
    }

    ~PqSpace() {}
};

}
