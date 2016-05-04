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
#include <vector>
#include <tuple>

///////////////////////////////////////////////////////////////////////////////
// host API
///////////////////////////////////////////////////////////////////////////////

#include "host.hpp"

///////////////////////////////////////////////////////////////////////////////
// lockstep measures: L_2 norm induced Euclidean metric
///////////////////////////////////////////////////////////////////////////////

double lockstepEuclidean1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);

    return lockstep(series0, length0,
                    series1, length1,
                    metric_euclidean_fixed<int, 1>());
}

double lockstepEuclidean2d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 2 == 0);

    return lockstep(series0, length0/2,
                    series1, length1/2,
                    metric_euclidean_fixed<int, 2>()) ;
}

double lockstepEuclidean3d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 3 == 0);

    return lockstep(series0, length0/3,
                    series1, length1/3,
                    metric_euclidean_fixed<int, 3>()) ;
}

double lockstepEuclidean4d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 4 == 0);

    return lockstep(series0, length0/4,
                    series1, length1/4,
                    metric_euclidean_fixed<int, 4>()) ;
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

    return lockstep(series0, length0/stride,
                    series1, length1/stride,
                    metric_euclidean_multivariate<int>(stride));
}

float lockstepEuclidean1f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);

    return lockstep(series0, length0,
                    series1, length1,
                    metric_euclidean_fixed<int, 1>()) ;
}

float lockstepEuclidean2f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 2 == 0);

    return lockstep(series0, length0/2,
                    series1, length1/2,
                    metric_euclidean_fixed<int, 2>()) ;
}

float lockstepEuclidean3f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 3 == 0);

    return lockstep(series0, length0/3,
                    series1, length1/3,
                    metric_euclidean_fixed<int, 3>()) ;
}

float lockstepEuclidean4f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 4 == 0);

    return lockstep(series0, length0/4,
                    series1, length1/4,
                    metric_euclidean_fixed<int, 4>()) ;
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

    return lockstep(series0, length0/stride,
                    series1, length1/stride,
                    metric_euclidean_multivariate<int>(stride));
}

///////////////////////////////////////////////////////////////////////////////
// lockstep measures: L_1 norm induced Manhattan metric
///////////////////////////////////////////////////////////////////////////////

double lockstepManhattan1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);

    return lockstep(series0, length0,
                    series1, length1,
                    metric_manhattan_fixed<int, 1>()) ;
}

double lockstepManhattan2d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 2 == 0);

    return lockstep(series0, length0/2,
                    series1, length1/2,
                    metric_manhattan_fixed<int, 2>()) ;
}

double lockstepManhattan3d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 3 == 0);

    return lockstep(series0, length0/3,
                    series1, length1/3,
                    metric_manhattan_fixed<int, 3>()) ;
}

double lockstepManhattan4d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 4 == 0);

    return lockstep(series0, length0/4,
                    series1, length1/4,
                    metric_manhattan_fixed<int, 4>()) ;
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

    return lockstep(series0, length0/stride,
                    series1, length1/stride,
                    metric_manhattan_multivariate<int>(stride));
}

float lockstepManhattan1f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);

    return lockstep(series0, length0,
                    series1, length1,
                    metric_manhattan_fixed<int, 1>()) ;
}

float lockstepManhattan2f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 2 == 0);

    return lockstep(series0, length0/2,
                    series1, length1/2,
                    metric_manhattan_fixed<int, 2>()) ;
}

float lockstepManhattan3f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 3 == 0);

    return lockstep(series0, length0/3,
                    series1, length1/3,
                    metric_manhattan_fixed<int, 3>()) ;
}

float lockstepManhattan4f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    // sanity checks
    assert(length0 == length1);
    assert(length0 % 4 == 0);

    return lockstep(series0, length0/4,
                    series1, length1/4,
                    metric_manhattan_fixed<int, 4>()) ;
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

    return lockstep(series0, length0/stride,
                    series1, length1/stride,
                    metric_manhattan_multivariate<int>(stride));
}

///////////////////////////////////////////////////////////////////////////////
// elastic measures: Euclidean-flavoured DTW similarity measure
///////////////////////////////////////////////////////////////////////////////

double elasticEuclideanDTW1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    std::vector<std::pair<int, int> > wpath;
    typedef metric_euclidean_fixed<int, 1> dist;
    return elastic_dtw<int, double, dist, 0, 0>(series0, length0,
                                                series1, length1,
                                                dist(), 0, wpath);
}

double elasticEuclideanCDTW1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      window) {

    std::vector<std::pair<int, int> > wpath;
    typedef metric_euclidean_fixed<int, 1> dist;
    return elastic_dtw<int, double, dist, 1, 0>(series0, length0,
                                                series1, length1,
                                                dist(), window, wpath);
}

double elasticEuclideanDTW1dBacktrace(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    std::vector<std::pair<int, int> > & wpath) {

    typedef metric_euclidean_fixed<int, 1> dist;
    return elastic_dtw<int, double, dist, 0, 1>(series0, length0,
                                                series1, length1,
                                                dist(), 0, wpath);
}

double elasticEuclideanCDTW1dBacktrace(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      window,
    std::vector<std::pair<int, int> > & wpath) {

    typedef metric_euclidean_fixed<int, 1> dist;
    return elastic_dtw<int, double, dist, 1, 1>(series0, length0,
                                                series1, length1,
                                                dist(), window, wpath);
}

double elasticEuclideanCDTW2dBacktrace(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      window,
    std::vector<std::pair<int, int> > & wpath) {

    // sanity checks
    assert(length0 % 2 == 0);
    assert(length1 % 2 == 0);

    typedef metric_euclidean_fixed<int, 2> dist;
    return elastic_dtw<int, double, dist, 1, 1>(series0, length0/2,
                                                series1, length1/2,
                                                dist(), window, wpath);
}
