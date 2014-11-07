#ifndef CPU_HELPER_HPP
#define CPU_HELPER_HPP

namespace helper {

    template <class index_t, class value_t>
    value_t lp_pow(value_t x, index_t p) {
        
        if (p < 0)
            return lp_pow(static_cast<value_t>(1)/x, -p);
        
        value_t result = static_cast<value_t>(1);
        
        while (p > 0) {
            
            if (p & 1) {
                result *= x;
            }
            
            x *= x;
            p >>= 1;
        }
        
        return result <  0 ? -result : result;
    }
}

#endif
