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

%module host
%{
/* Includes the header in the wrapper code */
#define SWIG_FILE_WITH_INIT
#include "host.hpp"
%}

%include "typemaps.i"
%include "numpy.i"
%include "std_vector.i"
%include "std_pair.i"

%init %{
import_array();
%}

namespace std {
    %template(WarpingNode) pair<int, int>;
    %template(WarpingPath) vector<pair<int, int> >;
}

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* series0, int length0)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* series1, int length1)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* env_lower, int len_lower)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* env_upper, int len_upper)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* matrix, int dimen0, int dimen1)};

%apply (float* INPLACE_ARRAY1, int DIM1) {(float* series0, int length0)};
%apply (float* INPLACE_ARRAY1, int DIM1) {(float* series1, int length1)};
%apply (float* INPLACE_ARRAY1, int DIM1) {(float* env_lower, int len_lower)};
%apply (float* INPLACE_ARRAY1, int DIM1) {(float* env_upper, int len_upper)};
%apply (float* INPLACE_ARRAY2, int DIM1, int DIM2) {(float* matrix, int dimen0, int dimen1)};

%include "host.hpp"
