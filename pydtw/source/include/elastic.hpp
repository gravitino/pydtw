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

#include <vector>

///////////////////////////////////////////////////////////////////////////////
// elastic measures
///////////////////////////////////////////////////////////////////////////////

template <
    typename index_t,
    typename value_t,
    typename funct_t>
value_t elastic_dtw_multivariate(
    value_t * series0,
    index_t   length0,
    value_t * series1,
    index_t   length1,
    funct_t   metric) {

    const index_t lane = length1+1;
    const index_t area = lane*(length0+1);
    std::vector<value_t> penalty(area, 0);
    

    return penalty[area-1];
}

template <
    typename index_t,
    typename value_t,
    typename funct_t>
value_t elstic_dtw_fixed(
    value_t * series0,
    index_t   length0,
    value_t * series1,
    index_t   length1,
    funct_t   metric) {

    value_t result = 0;

    return result;
}

#endif
