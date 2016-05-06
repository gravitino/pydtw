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

#ifndef PYDTW_RESIDUES_HPP
#define PYDTW_RESIDUES_HPP

///////////////////////////////////////////////////////////////////////////////
// optional OpenMP support
///////////////////////////////////////////////////////////////////////////////

#if defined(PYDTW_ENABLE_OPENMP)
#include <omp.h>
#endif

///////////////////////////////////////////////////////////////////////////////
// residues computation
///////////////////////////////////////////////////////////////////////////////

template <
    typename index_t,
    typename value_t,
    typename funct_t>
void residues_matrix(
    value_t * series0,
    value_t * series1,
    value_t * matrix,
    index_t   dimen0,
    index_t   dimen1,
    funct_t   metric) {

    #if defined(PYDTW_ENABLE_OPENMP)
    #pragma omp parallel for collapse(2)
    #endif
    for (index_t i = 0; i < dimen0; i++)
        for (index_t j = 0; j < dimen1; j++)
            matrix[i*dimen1+j] = metric(series0 + i*metric.stride,
                                        series1 + j*metric.stride);
}

#endif
