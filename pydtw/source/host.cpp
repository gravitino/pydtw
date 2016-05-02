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
// host API
//////////////////////////////////////////////////////////////////////////////

#include "host.hpp"

double lockstepEuclidean1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    return lockstep_fixed<int, 1>(series0, length0,
                                  series1, length1,
                                  metric_euclidean_fixed()) ;
}

float lockstepEuclidean1f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    return lockstep_fixed<int, 1>(series0, length0,
                                  series1, length1,
                                  metric_euclidean_fixed()) ;
}

double lockstepManhattan1d(
    double * series0,
    int      length0,
    double * series1,
    int      length1) {

    return lockstep_fixed<int, 1>(series0, length0,
                                  series1, length1,
                                  metric_manhattan_fixed()) ;
}

float lockstepManhattan1f(
    float * series0,
    int     length0,
    float * series1,
    int     length1) {

    return lockstep_fixed<int, 1>(series0, length0,
                                  series1, length1,
                                  metric_manhattan_fixed()) ;
}

/*
#include <algorithm>
#include <iostream>
#include <vector>

int main () {

    std::vector<float> series0(40, 0);
    std::vector<float> series1(40, 0);
    std::iota(series1.begin(), series1.end(), 1);

    auto metric_em = metric_euclidean_multivariate();
    std::cout << lockstep_multivariate(series0.data(), 10,
                                       series1.data(), 10,
                                       metric_em, 4) << std::endl;

    auto metric_mm = metric_manhattan_multivariate();
    std::cout << lockstep_multivariate(series0.data(), 10,
                                       series1.data(), 10,
                                       metric_mm, 4) << std::endl;

    auto metric_ef = metric_euclidean_fixed();
    std::cout << lockstep_fixed<int, 4>(series0.data(), 10,
                                        series1.data(), 10,
                                        metric_ef) << std::endl;

    auto metric_mf = metric_manhattan_fixed();
    std::cout << lockstep_fixed<int, 4>(series0.data(), 10,
                                        series1.data(), 10,
                                        metric_mf) << std::endl;

}
*/
