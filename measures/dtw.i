%module libdtw
%{
/* Includes the header in the wrapper code */
#include "dtw.hpp"
%}

%include "typemaps.i"
%include "std_vector.i"
%include "std_pair.i"

namespace std {
    %template(TimeSeries) vector<float>;
    %template(WarpingNode) pair<unsigned int, unsigned int>;
    %template(WarpingPath) vector<pair<unsigned int, unsigned int> >;
}

/* Parse the header file to generate wrappers */
%include "dtw.hpp"
