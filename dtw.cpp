#include "dtw.hpp"

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) < (y) ? (y) : (x))
#define min3(x, y, z) (min(x, min(y, z)))
#define sgn(x) ((x) < 0 ? -1 : 1)
#define abs(x) ((x)*sgn(x))

template <class itype, class ftype, bool squared> 
ftype meta_lp12(std::vector<ftype>* N, std::vector<ftype>* H) {

    // allover sum over residues
    ftype  sum = 0;
    
    for (itype i = 0; i < min(N->size(), H->size()); ++i) {
        ftype cost = N->at(i)-H->at(i);
            
        // Euclidean or Manhatten distance?
        if (squared)
            sum += cost*cost;
        else
            sum += abs(cost);
    }

    return sum;
}

template <class itype, class ftype, bool constrained, bool squared> 
ftype meta_dtw(std::vector<ftype>* N, std::vector<ftype>* H, itype w) {

    // vector sizes for convenience
    const itype Nsize = N->size(), Hsize = H->size();
    
    // alloc penalty matrix
    ftype *pen = new ftype[(Nsize+1)*(Hsize+1)];
  
    if (constrained) {
    
         // initialize fully penalty matrix
        for (itype i = 0; i < Nsize+1; ++i)
            for (itype j = 0; j < Hsize+1; ++j)
                pen[i*(Hsize+1)+j] = INFINITY;
        
        // adjust window if needed
        w = max(w, Nsize > Hsize ? Nsize-Hsize: Hsize-Nsize);
    
    } else {
  
        // initialize penalty matrix only in row and col zero 
        for (itype i = 0; i < Nsize+1; ++i)
            pen[i*(Hsize+1)] = INFINITY;
        for (itype j = 0; j < Hsize+1; ++j)
            pen[j] = INFINITY;
    }
    
    // set starting node to zero penalty
    pen[0] = 0;

    // relax penalty matrix
    for (itype i = 1; i < Nsize+1; ++i) {
    
        itype lower = 1, upper = Hsize;
    
        if (constrained) {
            lower = i > w ? i-w : 1;
            upper = min(Hsize, i+w);
        }
        
        for (itype j = lower; j < upper+1; ++j) {
            
            ftype cost = N->at(i-1)-H->at(j-1);
            
            // Euclidean or Manhatten distance?
            if (squared)
                cost *= cost;
            else
                cost *= sgn(cost);
            
            // backward relax three edges
            pen[i*(Hsize+1)+j] = cost+min3(pen[(i-1)*(Hsize+1)+j], 
                                           pen[(i-1)*(Hsize+1)+j-1], 
                                           pen[i*(Hsize+1)+j-1]);
        }
    }
    
    const ftype dist = pen[(Nsize+1)*(Hsize+1)-1];
    
    delete [] pen;

    return dist;
}

float euclidean(std::vector<float> *N, std::vector<float> *H) {
     
        return meta_lp12<unsigned int, float, true>(N, H);
}

float manhatten(std::vector<float> *N, std::vector<float> *H) {
     
        return meta_lp12<unsigned int, float, false>(N, H);
}

float dtw(std::vector<float> *N, std::vector<float> *H, bool squared) {
     
     if (squared)
        return meta_dtw<unsigned int, float, false, true>(N, H, 0);
     else
        return meta_dtw<unsigned int, float, false, false>(N, H, 0);
}

float cdtw(std::vector<float> *N, std::vector<float> *H, unsigned int w ,bool squared) {
     
     if (squared)
        return meta_dtw<unsigned int, float, true, true>(N, H, w);
     else
        return meta_dtw<unsigned int, float, true, false>(N, H, w);
}

