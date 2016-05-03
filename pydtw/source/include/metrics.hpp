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

#include <cmath>                      // std::min, std::acos
#include "qualifiers.hpp"             // qualifiers

///////////////////////////////////////////////////////////////////////////////
// local measures
///////////////////////////////////////////////////////////////////////////////

template <
    typename strde_t>
struct metric_euclidean_multivariate {

    const strde_t stride;
    metric_euclidean_multivariate(strde_t strde) : stride(strde) {}

    template <
        typename value_t> INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS
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

template <
    typename strde_t,
    strde_t  strde=1>
struct metric_euclidean_fixed {

    static constexpr strde_t stride = strde;

    template <
        typename value_t> INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS
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

template <
    typename strde_t>
struct metric_manhattan_multivariate {

    const strde_t stride;
    metric_manhattan_multivariate(strde_t strde) : stride(strde) {}

    template <
        typename value_t> INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS
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

template <
    typename strde_t,
    strde_t  strde=1>
struct metric_manhattan_fixed {

    static constexpr strde_t stride = strde;

    template <
        typename value_t> INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS
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

template <
    typename strde_t,
    strde_t  strde=4>
struct metric_quaternion_fixed {

    static constexpr strde_t stride = strde;

    template <
        typename value_t> INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS
    value_t operator()(
        const value_t * const __restrict__ entry0,
        const value_t * const __restrict__ entry1) const {

        value_t a_dot_b = 0;

        for (strde_t i = 0; i < stride; i++) {
            const value_t residue = entry0[i]*entry1[i];
            a_dot_b += residue;
        }
        a_dot_b = a_dot_b < 0 ? -a_dot_b : a_dot_b;
        return std::acos(std::min(value_t(1),a_dot_b));
    }
};

#endif
