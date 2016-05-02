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

#ifndef PYDTW_LOCKSTEP_HPP
#define PYDTW_LOCKSTEP_HPP

// TODO: maybe unroll (christian@metalabs.de)

template <
    typename strde_t,
    typename index_t,
    typename value_t,
    typename funct_t>
value_t lockstep_multivariate(
    value_t * series0,
    index_t   length0,
    value_t * series1,
    index_t   length1,
    funct_t   metric,
    strde_t   stride) {

    value_t result = 0;

    for (index_t i = 0; i < length0*stride; i += stride)
        result += metric(series0+i, series1+i, stride);

    return result;
}

// TODO: maybe unroll (christian@metalabs.de)

template <
    typename strde_t,
    strde_t stride=1,
    typename index_t,
    typename value_t,
    typename funct_t>
value_t lockstep_fixed(
    value_t * series0,
    index_t   length0,
    value_t * series1,
    index_t   length1,
    funct_t   metric) {

    value_t result = 0;

    for (index_t i = 0; i < length0*stride; i += stride)
        result += metric.template operator()<strde_t, stride>(series0+i, series1+i);

    return result;
}

#endif
