#ifndef CPU_GLOBAL_MEASURES_HPP
#define CPU_GLOBAL_MEASURES_HPP

#include <immintrin.h>   // avx instructions
#include "helper.hpp"    // convenience stuff

namespace lockstep {

    template <class index_t, class value_t>
    value_t dist_euclidean (value_t * query  , index_t M,
                            value_t * subject, index_t N) {
    
        assert (M == N);
        
        value_t result = static_cast<value_t>(0);
        
        for (index_t i = 0; i < N; ++i) {
            const value_t res = query[i] - subject[i];
            result += res*res;
        }
        
        return result;
    }
    
    template <class index_t, class value_t>
    value_t dist_manhattan (value_t * query  , index_t M,
                            value_t * subject, index_t N) {
    
        assert (M == N);
        
        value_t result = static_cast<value_t>(0);
        
        for (index_t i = 0; i < N; ++i) {
            const value_t res = query[i] - subject[i];
            result += res < 0 ? -res : res;
        }
        
        return result;
    }
    
    template <class index_t, class value_t> 
    value_t dist_lp (value_t * query  , index_t M,
                     value_t * subject, index_t N, index_t p) {
        
        assert (M == N);
        
        value_t result = static_cast<value_t>(0);
        
        for (index_t i = 0; i < N; ++i) {
            const value_t res = query[i] - subject[i];
            result += helper::lp_pow(res, p);
        }
        
        return result;
    }
    
    template <class index_t>
    double dist_euclidean_avx_d (double * query  , index_t M,
                                 double * subject, index_t N) {
    
        assert (M == N);
        
        double result = 0.0;
        const index_t border = (N/4)*4;
        __m256d cache = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
            
        for (index_t i = 0; i < border; i += 4) {
            __m256d a = _mm256_loadu_pd(query+i);
            __m256d b = _mm256_loadu_pd(subject+i);
            
            __m256d r = _mm256_sub_pd(a, b);
            cache = _mm256_add_pd(cache, _mm256_mul_pd(r, r));
        }
        
        __m256d hsum = _mm256_add_pd(cache, 
                       _mm256_permute2f128_pd(cache, cache, 0x1));
        _mm_store_sd(&result, _mm_hadd_pd(_mm256_castpd256_pd128(hsum),
                                          _mm256_castpd256_pd128(hsum)));
            
        for (index_t i = border; i < N; ++i) {
            const double res = query[i] - subject[i];
            result += res*res;
        }
        
        return result;
    }
    
    template <class index_t>
    float dist_euclidean_avx_f (float * query  , index_t M,
                                float * subject, index_t N) {
                                  
        assert (M == N);
        
        float result = 0.0;
        const index_t border = (N/8)*8;
        __m256 cache = _mm256_set_ps(0.0f, 0.0f, 0.0f, 0.0f,
                                     0.0f, 0.0f, 0.0f, 0.0f);
        
        for (index_t i = 0; i < border; i += 8) {
            __m256 r = _mm256_sub_ps(_mm256_loadu_ps(query+i),
                                     _mm256_loadu_ps(subject+i));
            cache = _mm256_add_ps(cache, _mm256_mul_ps(r, r));
        }

        __m256 hsum = _mm256_hadd_ps(cache, cache);
        hsum = _mm256_add_ps(hsum, _mm256_permute2f128_ps(hsum, hsum, 0x1));
        _mm_store_ss(&result, _mm_hadd_ps(_mm256_castps256_ps128(hsum), 
                                          _mm256_castps256_ps128(hsum)));
            
        for (index_t i = border; i < N; ++i) {
           const float res = query[i] - subject[i];
           result += res*res;
        }
        
        return result;
    }
    
    template <class index_t, class value_t>
    value_t dist_keogh_cid (value_t * query  , index_t M,
                            value_t * subject, index_t N) {
        
        assert (M == N);
        
        value_t CEQ = static_cast<value_t>(0),
                CES = static_cast<value_t>(0);
        
        for (index_t i = 0; i < M-1; ++i) {
            value_t res = query[i+1] - query[i];
            CEQ += res*res;
            
            res = subject[i+1] - subject[i];
            CES += res*res;
        }
        
        const value_t x = CEQ > CES ? CEQ : CES;
        const value_t y = CEQ < CES ? CEQ : CES;

        return dist_euclidean(query, M, subject, N) * ((x == y) ? 1 : x / y);
    }
}

#endif
