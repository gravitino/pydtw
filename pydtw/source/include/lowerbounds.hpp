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
#include "constants.hpp"
#include "qualifiers.hpp"
#include "mathematics.hpp"

///////////////////////////////////////////////////////////////////////////////
// streamed windowed min max envelope
///////////////////////////////////////////////////////////////////////////////

template <
    typename strde_t,
    typename index_t,
    typename value_t> INLINE_QUALIFIERS
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

template <
    typename index_t,
    typename value_t,
    typename funct_t> INLINE_QUALIFIERS ARCHITECTURE_QUALIFIERS
value_t elastic_LB_Kim4(
    value_t * series0,
    index_t   length0,
    value_t * series1,
    index_t   length1,
    funct_t   metric) {

    // stride
    const auto stride = metric.stride;

    //four entries of query and subject
    value_t * q0 = series0+0*stride;
    value_t * q1 = series0+1*stride;
    value_t * q2 = series0+2*stride;
    value_t * q3 = series0+3*stride;
    value_t * s0 = series1+0*stride;
    value_t * s1 = series1+1*stride;
    value_t * s2 = series1+2*stride;
    value_t * s3 = series1+3*stride;

    // linear memory penalty matrix
    value_t p00, p01, p02, p03,
            p10, p11, p12, p13, bsf, rst;

    // relax first row
    p00 = metric(q0, s0);
    p01 = p00 + metric(q0, s1);
    p02 = p01 + metric(q0, s2);
    p03 = p02 + metric(q0, s3);
    bsf = p03;

    // relax second row
    p10 = p00 + metric(q1, s0);
    p11 = pydtw_min3(p00, p01, p10) + metric(q1, s1);
    p12 = pydtw_min3(p01, p02, p11) + metric(q1, s2);
    p13 = pydtw_min3(p02, p03, p12) + metric(q1, s3);
    bsf = pydtw_min(bsf, p13);

    // relax third row
    p00 = p10 + metric(q2, s0);
    p01 = pydtw_min3(p10, p11, p00) + metric(q2, s1);
    p02 = pydtw_min3(p11, p12, p01) + metric(q2, s2);
    p03 = pydtw_min3(p12, p13, p02) + metric(q2, s3);
    bsf = pydtw_min(bsf, p03);

    // relax fourth row
    p10 = p00 + metric(q3, s0);
    p11 = pydtw_min3(p00, p01, p10) + metric(q3, s1);
    p12 = pydtw_min3(p01, p02, p11) + metric(q3, s2);
    p13 = pydtw_min3(p02, p03, p12) + metric(q3, s3);
    rst = pydtw_min(bsf, pydtw_min4(p10, p11, p12, p13));

    // four entries of query and subject
    q0 = series0+(length0-1)*stride;
    q1 = series0+(length0-2)*stride;
    q2 = series0+(length0-3)*stride;
    q3 = series0+(length0-4)*stride;
    s0 = series1+(length1-1)*stride;
    s1 = series1+(length1-2)*stride;
    s2 = series1+(length1-3)*stride;
    s3 = series1+(length1-4)*stride;

    // relax first row
    p00 = metric(q0, s0);
    p01 = p00 + metric(q0, s1);
    p02 = p01 + metric(q0, s2);
    p03 = p02 + metric(q0, s3);
    bsf = p03;

    // relax second row
    p10 = p00 + metric(q1, s0);
    p11 = pydtw_min3(p00, p01, p10) + metric(q1, s1);
    p12 = pydtw_min3(p01, p02, p11) + metric(q1, s2);
    p13 = pydtw_min3(p02, p03, p12) + metric(q1, s3);
    bsf = pydtw_min(bsf, p13);

    // relax third row
    p00 = p10 + metric(q2, s0);
    p01 = pydtw_min3(p10, p11, p00) + metric(q2, s1);
    p02 = pydtw_min3(p11, p12, p01) + metric(q2, s2);
    p03 = pydtw_min3(p12, p13, p02) + metric(q2, s3);
    bsf = pydtw_min(bsf, p03);

    // relax fourth row
    p10 = p00 + metric(q3, s0);
    p11 = pydtw_min3(p00, p01, p10) + metric(q3, s1);
    p12 = pydtw_min3(p01, p02, p11) + metric(q3, s2);
    p13 = pydtw_min3(p02, p03, p12) + metric(q3, s3);
    rst = pydtw_min(bsf, pydtw_min4(p10, p11, p12, p13)) + rst;

    return rst;
}


#endif
