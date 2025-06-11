#pragma once
#include "hnswlib.h"

namespace hnswlib {

template <typename Tdist, typename Tcorr>
static Tdist
InnerProductRef(const Tcorr* a, const Tcorr* b, size_t d) {
    Tdist sum = 0;
    for (size_t i = 0; i < d; i++) {
        sum += Tdist(a[i]) * Tdist(b[i]);
    }
    return sum;
}


int fvec_inner_product_int8_avx(const int8_t* x, const int8_t* y, size_t d) {
    __m256i msum256 = _mm256_setzero_si256();
    while (d >= 16) {
        __m256i ma = _mm256_cvtepi8_epi16(_mm_loadu_si128((const __m128i_u*)x));
        x += 16;
        __m256i mb = _mm256_cvtepi8_epi16(_mm_loadu_si128((const __m128i_u*)y));
        y += 16;
        msum256 = _mm256_add_epi32(msum256, _mm256_madd_epi16(ma, mb));
        d -= 16;
    }
    __m128i msum128 = _mm256_extracti128_si256(msum256, 1);
    msum128 = _mm_add_epi32(msum128, _mm256_extracti128_si256(msum256, 0));
    msum128 = _mm_hadd_epi32(msum128, msum128);
    msum128 = _mm_hadd_epi32(msum128, msum128);
    int sum = _mm_cvtsi128_si32(msum128);
    return d ? sum + InnerProductRef<int, int8_t>(x, y, d) : sum;
}

static float
InnerProductDistFuncIP(const void* a, const void* b, const void* d, float scale2) {
    size_t dim = *((size_t*)d);
    return (float)(fvec_inner_product_int8_avx((const int8_t*)a, (const int8_t*)b, dim)) / scale2;
}


static float
InnerProductDistFunc(const void* a, const void* b, const void* d, float scale2) {
    return  1.0 - InnerProductDistFuncIP(a, b, d, scale2);
}

class SpaceInt8 : public SpaceInterface<float> {
    DISTFUNC<float> fstdistfunc_;
    size_t data_size_;
    size_t dim_;

 public:
    SpaceInt8(size_t dim) {
        fstdistfunc_ = InnerProductDistFunc;
        dim_ = dim;
        data_size_ = dim * sizeof(int8_t);
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

    ~SpaceInt8() {}
};

}  // namespace hnswlib
