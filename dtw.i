%module libdtw
%{
/* Includes the header in the wrapper code */
#include "dtw.hpp"
%}

%include "typemaps.i"
%include "std_vector.i"

namespace std {
    %template(TimeSeries) vector<float>;
}

/* Parse the header file to generate wrappers */
%include "dtw.hpp"
