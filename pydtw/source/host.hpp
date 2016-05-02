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

#ifndef PYDTW_HOST_HPP
#define PYDTW_HOST_HPP

///////////////////////////////////////////////////////////////////////////////
// includes
///////////////////////////////////////////////////////////////////////////////

// metrics
#include "include/metrics.hpp"

// lock-step-measures
#include "include/lockstep.hpp"

///////////////////////////////////////////////////////////////////////////////
// exports
///////////////////////////////////////////////////////////////////////////////

double lockstepEuclidean1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

float lockstepEuclidean1f(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

double lockstepManhattan1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

float lockstepManhattan1f(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

#endif
