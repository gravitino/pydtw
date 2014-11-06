#ifndef CPU_GLOBAL_MEASURES_HPP
#define CPU_GLOBAL_MEASURES_HPP

#include <iostream>
#include <type_traits>   // static type control
#include <immintrin.h>   // avx instructions

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
    value_t dist_euclidean_avx (value_t * query  , index_t M,
                                value_t * subject, index_t N) {
    
        assert (M == N);
        
        value_t result = static_cast<value_t>(0);
        
        if (std::is_same<float, value_t>()) {
        
            index_t border = (N/8)*8;
            __m256 cache;
        
            for (index_t i = 0; i < border; i += 8) {            
                __m256 r = _mm256_sub_ps(_mm256_loadu_ps(query+i),
                                         _mm256_loadu_ps(subject+i));
                cache = _mm256_add_ps(cache, _mm256_mul_ps(r, r));                
            }
            
            // I should rethink the next five lines
            const __m128 x128 = _mm_add_ps(_mm256_extractf128_ps(cache, 1), 
                                           _mm256_castps256_ps128(cache));
            const __m128 x64 = _mm_add_ps(x128, _mm_movehl_ps(x128, x128));
            const __m128 x32 = _mm_add_ss(x64, _mm_shuffle_ps(x64, x64, 0x55));
            result += _mm_cvtss_f32(x32);
            
            for (index_t i = border; i < N; ++i) {
                const value_t res = query[i] - subject[i];
                result += res*res;
            }
            
        } else if (std::is_same<double, value_t>()) {
            return -1;
        } else {
            return -1;
        }
        
        return result;
    }

}

#endif
