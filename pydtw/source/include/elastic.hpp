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

#ifndef PYDTW_ELASTIC_HPP
#define PYDTW_ELASTIC_HPP

///////////////////////////////////////////////////////////////////////////////
// optional OpenMP support
///////////////////////////////////////////////////////////////////////////////

#if defined(PYDTW_ENABLE_OPENMP)
#include <omp.h>
#endif

///////////////////////////////////////////////////////////////////////////////
// includes
///////////////////////////////////////////////////////////////////////////////

#include "constants.hpp"

///////////////////////////////////////////////////////////////////////////////
// elastic measures
///////////////////////////////////////////////////////////////////////////////

template <
    typename index_t,
    typename value_t,
    typename funct_t,
    bool compute_backtrace=false>
value_t elastic_dtw(
    value_t * series0,
    index_t   length0,
    value_t * series1,
    index_t   length1,
    funct_t   metric) {

    // convenience variables
    const index_t lane_i = length0+1;
    const index_t lane_j = length1+1;
    const index_t area = lane_i*lane_j;

    value_t * penalty = new value_t[area];
    uint8_t * predecs = nullptr;

    // backtracing needs quadratic memory: suppress if not needed
    if (compute_backtrace)
        predecs = new uint8_t[area];

    // initialize penalty matrix
    for (index_t j = 1; j < lane_j; j++)
        penalty[0*lane_j+j] = PYDTW_CONSTANTS_INFINITY;
    for (index_t i = 1; i < lane_i; i++)
        penalty[i*lane_j+0] = PYDTW_CONSTANTS_INFINITY;
    penalty[0] = 0;

    // start relaxing
    for (index_t i = 1; i < lane_i; i++) {

        // compute bounds for the window
        const index_t lower = 1;
        const index_t upper = lane_j;

        for (index_t j = lower; j < upper; j++) {

            // compute local measure
            const value_t cost = metric(series0+(i-1)*metric.stride,
                                        series1+(j-1)*metric.stride);

            // cache matrix entries
            const value_t diag = penalty[(i-1)*lane_j+(j-1)];
            const value_t abve = penalty[(i-1)*lane_j+(j+0)];
            const value_t left = penalty[(i-0)*lane_j+(j-1)];

            // prefer diagonal steps if tied
            value_t bsf_value = diag;
            index_t bsf_index = PYDTW_CONSTANTS_DIAGONAL;

            // search best extension of the warping path
            if (bsf_value > abve) {
                bsf_value = abve;
                if (compute_backtrace)
                    bsf_index = PYDTW_CONSTANTS_ABOVE;
            }
            if (bsf_value > left) {
                bsf_value = left;
                if (compute_backtrace)
                    bsf_index = PYDTW_CONSTANTS_LEFT;
            }

            // relax cell
            penalty[i*lane_j+j] = cost + bsf_value;
        }
    }

    const value_t result = penalty[area-1];
    delete [] penalty;

    if (compute_backtrace)
        delete [] predecs;

    return result;
}

#endif
