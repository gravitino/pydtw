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

#include <algorithm>                     // std::reverse
#include <cstdint>                       // uint8_t
#include <vector>                        // std::vector
#include "constants.hpp"                 // constants
#include "lowerbounds.hpp"               // lemire_min_max 

///////////////////////////////////////////////////////////////////////////////
// elastic measures
///////////////////////////////////////////////////////////////////////////////

template <
    typename index_t,
    typename value_t,
    typename funct_t,
    bool sakoe_constrained=false,
    bool compute_backtrace=false>
value_t elastic_dtw(
    value_t * series0,
    index_t   length0,
    value_t * series1,
    index_t   length1,
    funct_t   metric,
    index_t   window,
    std::vector<std::pair<index_t, index_t> >& wpath) {

    // Sakoe-Chiba-constrained DTW only for same length
    if (sakoe_constrained)
        assert(length0 == length1);

    // convenience variables
    const index_t lane_i = length0+1;
    const index_t lane_j = length1+1;
    const index_t area = lane_i*lane_j;

    // penalty and predecessor matrix
    value_t * penalty = new value_t[2*lane_j];
    uint8_t * predecs = nullptr;

    // initialize penalty matrix in linear memory: that's fine
    for (index_t j = 1; j < lane_j; j++)
        penalty[0*lane_j+j] = PYDTW_CONSTANTS_INFINITY;
    penalty[0] = 0;

    // backtracing needs quadratic memory: suppress if not needed
    if (compute_backtrace) {
        predecs = new uint8_t[area];
        for (index_t i = 1; i < lane_i; i++)
            predecs[i*lane_j+0] = PYDTW_CONSTANTS_NO_SOURCE;
        for (index_t j = 0; j < lane_j; j++)
            predecs[0*lane_j+j] = PYDTW_CONSTANTS_NO_SOURCE;
    }

    // start relaxing
    for (index_t i = 1; i < lane_i; i++) {

        // linear memory indexing
        const index_t this_lane = i & 1;
        const index_t prev_lane = !this_lane;

        // compute bounds for the window
        index_t lower, upper;
        if (sakoe_constrained) {
            lower = pydtw_max(i-window, 1);
            upper = pydtw_min(i+window+1, lane_j);
        } else {
            lower = 1;
            upper = lane_j;
        }

        // set left cell to INFINITY
        penalty[this_lane*lane_j+lower-1] = PYDTW_CONSTANTS_INFINITY;

        for (index_t j = lower; j < upper; j++) {

            // compute local measure
            const value_t cost = metric(series0+(i-1)*metric.stride,
                                        series1+(j-1)*metric.stride);

            // cache matrix entries
            const value_t diag = penalty[prev_lane*lane_j+(j-1)];
            const value_t abve = penalty[prev_lane*lane_j+(j+0)];
            const value_t left = penalty[this_lane*lane_j+(j-1)];

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
            penalty[this_lane*lane_j+j] = cost + bsf_value;

            // remember predecessor
            if (compute_backtrace)
                predecs[i*lane_j+j] = bsf_index;
        }

        // set right cell to INFINITY
        if (upper < lane_j)
            penalty[this_lane*lane_j+upper] = PYDTW_CONSTANTS_INFINITY;
    }

    // get the last lane, fetch the last cell, and free penalty matrix
    const index_t last_lane = !(lane_i & 1);
    const value_t result = penalty[last_lane*lane_j+lane_j-1];
    delete [] penalty;

    // if backtracing enabled follow path
    if (compute_backtrace) {

        index_t i = lane_i-1, j = lane_j-1;
        uint8_t direction = predecs[i*lane_j+j];
        wpath.push_back(std::pair<index_t, index_t>(i-1, j-1));

        while (direction != PYDTW_CONSTANTS_NO_SOURCE) {

            if (direction == PYDTW_CONSTANTS_DIAGONAL) {
                i -= 1;
                j -= 1;
            }

            if (direction == PYDTW_CONSTANTS_ABOVE) {
                i -= 1;
            }

            if (direction == PYDTW_CONSTANTS_LEFT) {
                j -= 1;
            }

            // update predecessor information and warping path
            direction = predecs[i*lane_j+j];
            wpath.push_back(std::pair<index_t, index_t>(i-1, j-1));
        }

        // remove the last added node, reverse result and free memory
        wpath.pop_back();
        std::reverse(wpath.begin(), wpath.end());
        delete [] predecs;
    }

    return result;
}

#endif
