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

#ifndef PYDTW_LOWERBOUNDS_HPP
#define PYDTW_LOWERBOUNDS_HPP

///////////////////////////////////////////////////////////////////////////////
// includes
///////////////////////////////////////////////////////////////////////////////

#include <deque>

///////////////////////////////////////////////////////////////////////////////
// streamed windowed min max envelope
///////////////////////////////////////////////////////////////////////////////

template <
    typename strde_t,
    typename index_t,
    typename value_t>
void elastic_lemire_min_max(
    value_t * series0,
    value_t * env_lower,
    value_t * env_upper,
    index_t   length,
    index_t   window,
    strde_t   stride) {

    // this is note cache optimal unless we perform an AOS2SOA conversion
    for (index_t dim = 0; dim < stride; dim++) {

        std::deque<index_t> indices_lower;
        std::deque<index_t> indices_upper;

        indices_lower.push_back(0);
        indices_upper.push_back(0);

        for (index_t i = 1; i < length; i++) {

            const index_t index = (i-window-1)*stride+dim;
            if (i > window) {
                env_upper[index] = series0[indices_upper.front()*stride+dim];
                env_lower[index] = series0[indices_lower.front()*stride+dim];
            }

            const value_t ying = series0[(i-0)*stride+dim];
            const value_t yang = series0[(i-1)*stride+dim];
            if (ying > yang) {
                indices_upper.pop_back();
                while (!indices_upper.empty() &&
                       ying > series0[indices_upper.back()*stride+dim])
                    indices_upper.pop_back();
            } else {
                indices_lower.pop_back();
                while (!indices_lower.empty() &&
                       ying < series0[indices_lower.back()*stride+dim])
                    indices_lower.pop_back();
            }

            indices_lower.push_back(i);
            indices_upper.push_back(i);

            const index_t base = 2*window+1;
            if (i == base+indices_upper.front())
                indices_upper.pop_front();
            else if (i == base+indices_lower.front())
                indices_lower.pop_front();
        }

        for (index_t i = length; i < length+window+1; i++) {

            const index_t index = (i-window-1)*stride+dim;
            env_lower[index] = series0[indices_lower.front()*stride+dim];
            env_upper[index] = series0[indices_upper.front()*stride+dim];

            const index_t base = 2*window+1;
            if (i-indices_upper.front() >= base)
                indices_upper.pop_front();
            if (i-indices_lower.front() >= base)
                indices_lower.pop_front();
        }
    }
}

#endif
