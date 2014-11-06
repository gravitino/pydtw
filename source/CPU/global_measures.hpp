#ifndef CPU_GLOBAL_MEASURES_HPP
#define CPU_GLOBAL_MEASURES_HPP

#include <iostream>

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

}

#endif
