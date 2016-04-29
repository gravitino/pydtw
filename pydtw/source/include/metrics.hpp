///////////////////////////////////////////////////////////////////////////////
//    This file is part of pydtw.
//
//    pydtw is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    pydtw is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with pydtw.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////////

#ifndef PYDTW_METRICS_HPP
#define PYDTW_METRICS_HPP

///////////////////////////////////////////////////////////////////////////////
// includes
///////////////////////////////////////////////////////////////////////////////

#include <include/qualifiers.hpp>             // qualifiers

///////////////////////////////////////////////////////////////////////////////
// local measures
///////////////////////////////////////////////////////////////////////////////

struct metric_euclidean_multivariate {
    template <
        typename strde_t,
        typename value_t> 

    INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS

    value_t operator()(
        const value_t * const __restrict__ entry0, 
        const value_t * const __restrict__ entry1,
        const strde_t stride=1) const {

        value_t result = 0;
        
        for (strde_t i = 0; i < stride; i++) {
            const value_t residue = entry0[i]-entry1[i];
            result += residue*residue;
        }

        return result;
    }
};

struct metric_euclidean_fixed {

    template <
        typename strde_t,
        strde_t stride=1,
        typename value_t> 

    INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS

    value_t operator()(
        const value_t * const __restrict__ entry0, 
        const value_t * const __restrict__ entry1) const {

        value_t result = 0;
        
        for (strde_t i = 0; i < stride; i++) {
            const value_t residue = entry0[i]-entry1[i];
            result += residue*residue;
        }

        return result;
    }
};

struct metric_manhattan_multivariate {

    template <
        typename strde_t,
        typename value_t> 

    INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS

    value_t operator()(
        const value_t * const __restrict__ entry0, 
        const value_t * const __restrict__ entry1,
        const strde_t stride=1) const {

        value_t result = 0;
        
        for (strde_t i = 0; i < stride; i++) {
            const value_t residue = entry0[i]-entry1[i];
            result += residue < 0 ? -residue : residue;
        }

        return result;
    }
};

struct metric_manhattan_fixed {

    template <
        typename strde_t,
        strde_t stride=1,
        typename value_t> 

    INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS

    value_t operator()(
        const value_t * const __restrict__ entry0, 
        const value_t * const __restrict__ entry1) const {

        value_t result = 0;
        
        for (strde_t i = 0; i < stride; i++) {
            const value_t residue = entry0[i]-entry1[i];
            result += residue < 0 ? -residue : residue;
        }

        return result;
    }
};

#endif
