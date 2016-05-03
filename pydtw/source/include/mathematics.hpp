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

#ifndef PYDTW_MATHEMATICS_HPP
#define PYDTW_MATHEMATICS_HPP

///////////////////////////////////////////////////////////////////////////////
// includes
///////////////////////////////////////////////////////////////////////////////

#include "qualifiers.hpp"

template <
    typename value_t> INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS
value_t pydtw_abs(const value_t& x) {
    return x < 0 ? -x : x;
}

template <
    typename value_t> INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS
value_t pydtw_min(const value_t& x, const value_t& y) {
    return x < y ? x : y;
}

template <
    typename value_t> INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS
value_t pydtw_max(const value_t& x, const value_t& y) {
    return x > y ? x : y;
}


#endif
