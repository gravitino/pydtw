%module libglobalmeasures

%{
/* Includes the header in the wrapper code */
#define SWIG_FILE_WITH_INIT
#include "global_measures.hpp"
%}

%include "typemaps.i"
%include "numpy.i"

%init %{
import_array();
%}

/* Parse the header file to generate wrappers */

%include "global_measures.hpp"

// DOUBLE PRECISION
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* query,   int M)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* subject, int N)};
%template(dist_euclidean_d) lockstep::dist_euclidean<int, double>;
%template(dist_manhattan_d) lockstep::dist_manhattan<int, double>;

// SINGLE PRECISION
%apply (float* INPLACE_ARRAY1, int DIM1) {(float* query,   int M)};
%apply (float* INPLACE_ARRAY1, int DIM1) {(float* subject, int N)};
%template(dist_euclidean_f) lockstep::dist_euclidean<int, float>;
%template(dist_manhattan_f) lockstep::dist_manhattan<int, float>;

%template(dist_euclidean_avx_f) lockstep::dist_euclidean_avx<int, float>;


