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

///////////////////////////////////////////////////////////////////////////////
// includes
///////////////////////////////////////////////////////////////////////////////

#include <assert.h>

///////////////////////////////////////////////////////////////////////////////
// host API
///////////////////////////////////////////////////////////////////////////////

#include "host.hpp"

///////////////////////////////////////////////////////////////////////////////
// lockstep measures: Euclidean
///////////////////////////////////////////////////////////////////////////////

double lockstepEuclidean1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);

    return lockstep_fixed<int, 1>(series0, length0,
                                  series1, length1,
                                  metric_euclidean_fixed()) ;
}

double lockstepEuclidean2d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 2 == 0);

    return lockstep_fixed<int, 2>(series0, length0/2,
                                  series1, length1/2,
                                  metric_euclidean_fixed()) ;
}

double lockstepEuclidean3d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 3 == 0);

    return lockstep_fixed<int, 3>(series0, length0/3,
                                  series1, length1/3,
                                  metric_euclidean_fixed()) ;
}

double lockstepEuclidean4d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 4 == 0);

    return lockstep_fixed<int, 4>(series0, length0/4,
                                  series1, length1/4,
                                  metric_euclidean_fixed()) ;
}

double lockstepEuclideanNd(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      stride) {

    // sanity checks
    assert(length0 == length1);
    assert(stride > 0);
    assert(length0 % stride == 0);

    return lockstep_multivariate(series0, length0/stride,
                                 series1, length1/stride,
                                 metric_euclidean_multivariate(), stride);
}

float lockstepEuclidean1f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);

    return lockstep_fixed<int, 1>(series0, length0,
                                  series1, length1,
                                  metric_euclidean_fixed()) ;
}

float lockstepEuclidean2f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 2 == 0);

    return lockstep_fixed<int, 2>(series0, length0/2,
                                  series1, length1/2,
                                  metric_euclidean_fixed()) ;
}

float lockstepEuclidean3f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 3 == 0);

    return lockstep_fixed<int, 3>(series0, length0/3,
                                  series1, length1/3,
                                  metric_euclidean_fixed()) ;
}

float lockstepEuclidean4f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 4 == 0);

    return lockstep_fixed<int, 4>(series0, length0/4,
                                  series1, length1/4,
                                  metric_euclidean_fixed()) ;
}

float lockstepEuclideanNf(
    float * series0,
    int     length0,
    float * series1,
    int     length1,
    int     stride) {

    // sanity checks
    assert(length0 == length1);
    assert(stride > 0);
    assert(length0 % stride == 0);

    return lockstep_multivariate(series0, length0/stride,
                                 series1, length1/stride,
                                 metric_euclidean_multivariate(), stride);
}

///////////////////////////////////////////////////////////////////////////////
// lockstep measures: Manhattan
///////////////////////////////////////////////////////////////////////////////

double lockstepManhattan1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);

    return lockstep_fixed<int, 1>(series0, length0,
                                  series1, length1,
                                  metric_manhattan_fixed()) ;
}

double lockstepManhattan2d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 2 == 0);

    return lockstep_fixed<int, 2>(series0, length0/2,
                                  series1, length1/2,
                                  metric_manhattan_fixed()) ;
}

double lockstepManhattan3d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 3 == 0);

    return lockstep_fixed<int, 3>(series0, length0/3,
                                  series1, length1/3,
                                  metric_manhattan_fixed()) ;
}

double lockstepManhattan4d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 4 == 0);

    return lockstep_fixed<int, 4>(series0, length0/4,
                                  series1, length1/4,
                                  metric_manhattan_fixed()) ;
}

double lockstepManhattanNd(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      stride) {

    // sanity checks
    assert(length0 == length1);
    assert(stride > 0);
    assert(length0 % stride == 0);

    return lockstep_multivariate(series0, length0/stride,
                                 series1, length1/stride,
                                 metric_manhattan_multivariate(), stride);
}


float lockstepManhattan1f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);

    return lockstep_fixed<int, 1>(series0, length0,
                                  series1, length1,
                                  metric_manhattan_fixed()) ;
}

float lockstepManhattan2f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 2 == 0);

    return lockstep_fixed<int, 2>(series0, length0/2,
                                  series1, length1/2,
                                  metric_manhattan_fixed()) ;
}

float lockstepManhattan3f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 3 == 0);

    return lockstep_fixed<int, 3>(series0, length0/3,
                                  series1, length1/3,
                                  metric_manhattan_fixed()) ;
}

float lockstepManhattan4f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 4 == 0);

    return lockstep_fixed<int, 4>(series0, length0/4,
                                  series1, length1/4,
                                  metric_manhattan_fixed()) ;
}

float lockstepManhattanNf(
    float * series0,
    int     length0,
    float * series1,
    int     length1,
    int     stride) {

    // sanity checks
    assert(length0 == length1);
    assert(stride > 0);
    assert(length0 % stride == 0);

    return lockstep_multivariate(series0, length0/stride,
                                 series1, length1/stride,
                                 metric_manhattan_multivariate(), stride);
}
