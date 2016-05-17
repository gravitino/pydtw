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

// lock-step measures
#include "include/lockstep.hpp"

// elastic measures
#include "include/elastic.hpp"

// residue matrix
#include "include/residues.hpp"

///////////////////////////////////////////////////////////////////////////////
// lockstep measures: L_2 norm induced Euclidean metric
///////////////////////////////////////////////////////////////////////////////

double lockstepEuclidean1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

double lockstepEuclidean2d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

double lockstepEuclidean3d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

double lockstepEuclidean4d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

double lockstepEuclideanNd(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      stride);

float lockstepEuclidean1f(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

float lockstepEuclidean2f(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

float lockstepEuclidean3f(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

float lockstepEuclidean4f(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

float lockstepEuclideanNf(
    float * series0,
    int     length0,
    float * series1,
    int     length1,
    int     stride);

///////////////////////////////////////////////////////////////////////////////
// lockstep measures: L_1 norm induced Manhattan metric
///////////////////////////////////////////////////////////////////////////////

double lockstepManhattan1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

double lockstepManhattan2d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

double lockstepManhattan3d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

double lockstepManhattan4d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

double lockstepManhattanNd(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      stride);

float lockstepManhattan1f(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

float lockstepManhattan2f(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

float lockstepManhattan3f(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

float lockstepManhattan4f(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

float lockstepManhattanNf(
    float * series0,
    int     length0,
    float * series1,
    int     length1,
    int     stride);

///////////////////////////////////////////////////////////////////////////////
// lockstep measures: quaternion metric
///////////////////////////////////////////////////////////////////////////////

float lockstepQuaternionf(
    float * series0,
    int     length0,
    float * series1,
    int     length1);

double lockstepQuaterniond(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

///////////////////////////////////////////////////////////////////////////////
// elastic measures: Euclidean-flavoured DTW similarity measure
///////////////////////////////////////////////////////////////////////////////

double elasticEuclideanDTW1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

double elasticEuclideanCDTW1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      window);

double elasticEuclideanDTW1dBacktrace(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    std::vector<std::pair<int, int> > & wpath);

double elasticEuclideanCDTW1dBacktrace(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      window,
    std::vector<std::pair<int, int> > & wpath);

double elasticEuclideanCDTW2dBacktrace(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      window,
    std::vector<std::pair<int, int> > & wpath);

///////////////////////////////////////////////////////////////////////////////
// elastic measures: quaternion DTW similarity measure
///////////////////////////////////////////////////////////////////////////////

double elasticQuaternionCDTWd(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    int      window);

///////////////////////////////////////////////////////////////////////////////
// elastic measures: LB_Keogh envelopes using Lemire efficient streamed min/max
///////////////////////////////////////////////////////////////////////////////

void elasticWarpingEnvelopeNd(
    double * series0,
    int      length0,
    double * env_lower,
    int      len_lower,
    double * env_upper,
    int      len_upper,
    int      window,
    int      stride);

///////////////////////////////////////////////////////////////////////////////
// elastic measures: LB_Kim
///////////////////////////////////////////////////////////////////////////////

double elasticLBKim4d(
    double * series0,
    int      length0,
    double * series1,
    int      length1);

///////////////////////////////////////////////////////////////////////////////
// elastic measures: LB_Keogh
///////////////////////////////////////////////////////////////////////////////

double elasticLBKeogh4d(
    double * env_lower,
    int      len_lower,
    double * env_upper,
    int      len_upper,
    double * series0,
    int      length0);

///////////////////////////////////////////////////////////////////////////////
// residue matrix for arbitrary local metrics
///////////////////////////////////////////////////////////////////////////////

void residuesMatrixEuclideanNd(
    double * series0,
    int      length0,
    double * series1,
    int      length1,
    double * matrix,
    int      dimen0,
    int      dimen1,
    int      stride);
    
#endif
