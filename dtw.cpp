#include "dtw.hpp"

#define min(x, y) (x < y ? x : y)
#define max(x, y) (x < y ? y : x)
#define min3(x, y, z) (min(x, min(y, z)))
#define sgn(x) (x < 0 ? -1 : 1)

template <class itype, class ftype, bool squared> 
ftype meta_dtw(std::vector<ftype>* N, std::vector<ftype>* H) {

    // vector sizes for convenience
    const itype Nsize = N->size(), Hsize = H->size();

    // alloc penalty matrix
    ftype *pen = new ftype[(Nsize+1)*(Hsize+1)];
  
    // initialize penalty matrix
    for (itype i = 1; i < Nsize+1; ++i)
        pen[i*(Hsize+1)] = INFINITY;
    for (itype j = 1; j < Hsize+1; ++j)
        pen[j] = INFINITY;
    pen[0] = 0;

    // relax penalty matrix
    for (itype i = 1; i < Nsize+1; ++i)
        for (itype j = 1; j < Hsize+1; ++j) {
            ftype cost = N->at(i-1)-H->at(j-1);
            
            // Euclidian or Manhatten distance?
            if (squared)
                cost *= cost;
            else
                cost *= sgn(cost);
            
            // backward relax three edges
            pen[i*(Hsize+1)+j] = cost+min3(pen[(i-1)*(Hsize+1)+j], 
                                           pen[(i-1)*(Hsize+1)+j-1], 
                                           pen[i*(Hsize+1)+j-1]);
        }

    const ftype dist = pen[(Nsize+1)*(Hsize+1)-1];

    // free memory
    delete [] pen;

    return dist;
}

float dtw(std::vector<float> *N, std::vector<float> *H, bool squared) {
     
     if (squared)
        return meta_dtw<unsigned int, float, true>(N, H);
     else
        return meta_dtw<unsigned int, float, false>(N, H);
}
