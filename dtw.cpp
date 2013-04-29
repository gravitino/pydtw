#include "dtw.hpp"

///////////////////////////////////////////////////////////////////////////////
// distance measures
///////////////////////////////////////////////////////////////////////////////

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

template <class itype, class ftype, bool constrained, bool squared, bool backtrace> 
ftype meta_dtw(std::vector<ftype>* N, std::vector<ftype>* H, itype w, 
               std::vector<std::pair<itype, itype> >* path) {

    // vector sizes for convenience
    const itype Nsize = N->size(), Hsize = H->size();
    
    // alloc penalty matrix
    ftype *pen = new ftype[(Nsize+1)*(Hsize+1)];
    char *pre = 0;
  
    if (backtrace) {
        pre = new char[(Nsize+1)*(Hsize+1)];
    } else {
        unused(pre);
    }
  
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
            
            if (backtrace) {
                // remember predecessor node
                pre[i*(Hsize+1)+j] = argmin3(pen[(i-1)*(Hsize+1)+j], 
                                             pen[(i-1)*(Hsize+1)+j-1], 
                                             pen[i*(Hsize+1)+j-1]);
            }
        }
    }
    
    // remember the optimal distance measure
    const ftype dist = pen[(Nsize+1)*(Hsize+1)-1];
    
    // free memory
    delete [] pen;
    
    if (backtrace) {
        
        // overwrite existing warping path
        path->clear();
        
        // start at the end of the warping path
        itype i = Nsize, j = Hsize;
        
        // until starting node reached
        while (i != 0 && j !=0) {
        
            // add node to warping path
            path->push_back(std::pair<itype, itype>(i-1, j-1));
            
            // manipulate indices of needle
            if (pre[i*(Hsize+1)+j] == 0) {
                i -= 1;
            }
            if (pre[i*(Hsize+1)+j] == 1) {
                i -= 1; j -=1;
            }
            if (pre[i*(Hsize+1)+j] == 2) {
                j -=1;
            }
        }

        // reverse the warping path
        std::reverse(path->begin(), path->end());
    
        // free memory
        delete [] pre;
    }
    
    return dist;
}

///////////////////////////////////////////////////////////////////////////////
// wrappers for distance measures
///////////////////////////////////////////////////////////////////////////////

float dist_euclidean(std::vector<float> *N, std::vector<float> *H) {
     
        return meta_lp12<unsigned int, float, true>(N, H);
}

float dist_manhatten(std::vector<float> *N, std::vector<float> *H) {
     
        return meta_lp12<unsigned int, float, false>(N, H);
}

float dist_dtw(std::vector<float> *N, std::vector<float> *H, bool squared) {
     
     if (squared)
        return meta_dtw<unsigned int, float, false, true, false>(N, H, 0, 0);
     else
        return meta_dtw<unsigned int, float, false, false, false>(N, H, 0, 0);
}

float dist_dtw_backtrace(std::vector<float> *N, std::vector<float> *H, 
      std::vector<std::pair<unsigned int, unsigned int> >* path, bool squared) {
     
     if (squared)
        return meta_dtw<unsigned int, float, false, true, true>(N, H, 0, path);
     else
        return meta_dtw<unsigned int, float, false, false, true>(N, H, 0, path);
}

float dist_cdtw(std::vector<float> *N, std::vector<float> *H, unsigned int w ,bool squared) {
     
     if (squared)
        return meta_dtw<unsigned int, float, true, true, false>(N, H, w, 0);
     else
        return meta_dtw<unsigned int, float, true, false, false>(N, H, w, 0);
}

float dist_cdtw_backtrace(std::vector<float> *N, std::vector<float> *H, unsigned int w, 
      std::vector<std::pair<unsigned int, unsigned int> >* path, bool squared) {
     
     if (squared)
        return meta_dtw<unsigned int, float, true, true, true>(N, H, w, path);
     else
        return meta_dtw<unsigned int, float, true, false, true>(N, H, w, path);
}

///////////////////////////////////////////////////////////////////////////////
// lower bound - related methods
///////////////////////////////////////////////////////////////////////////////

template <class itype, class ftype> 
int meta_envelope(std::vector<ftype> *series, itype w, 
                  std::vector<ftype> *L, std::vector<ftype> *U) {
    
    // remember n for convenience
    itype n = series->size();
    
    // envelope preinitialized with zeros
    L->clear();
    U->clear();
    L->resize(n, 0);
    U->resize(n, 0);
    
    // Daniel Lemire's windowed min-max algorithm in O(3n):
    std::list<itype> u = std::list<itype>();
    std::list<itype> l = std::list<itype>();

    u.push_back(0);
    l.push_back(0);

    for (itype i = 1; i < n; ++i) {

        if (i > w) {

            U->at(i-w-1) = series->at(u.front());
            L->at(i-w-1) = series->at(l.front());
        }
        
        if (series->at(i) > series->at(i-1)) {
            
            u.pop_back();
            while (!u.empty() && series->at(i) > series->at(u.back()))
                u.pop_back();
        } else {

            l.pop_back();
            while (!l.empty() && series->at(i) < series->at(l.back()))
                l.pop_back();
        }
        
        u.push_back(i);
        l.push_back(i);
        
        if (i == 2*w+1+u.front())
            u.pop_front();
        else if (i == 2*w+1+l.front())
            l.pop_front();
    }

    for (itype i = n; i < n+w+1; ++i) {

        U->at(i-w-1) = series->at(u.front());
        L->at(i-w-1) = series->at(l.front());

        if (i-u.front() >= 2*w+1)
            u.pop_front();

        if (i-l.front() >= 2*w+1)
            l.pop_front();
    }

    return 0;
}

template <class itype, class ftype, bool precalculated, bool onquery, bool squared> 
ftype meta_lb_keogh(std::vector<ftype> *query, std::vector<ftype> *subject, 
                    itype w, std::vector<ftype> *L, std::vector<ftype> *U) {

    ftype penalty = 0;

    // if envelope not calculated do it
    if (!precalculated) {
    
        L = new std::vector<ftype>();
        U = new std::vector<ftype>();
    
        if (onquery)
            meta_envelope<itype, ftype> (query, w, L, U);
        else
            meta_envelope<itype, ftype> (subject, w, L, U);
    } else {
        if (onquery)
            query = new std::vector<ftype>(subject->size(), 0);
        else
            subject = new std::vector<ftype>(query->size(), 0);
    }
    
    // calculate sum of differences to the given envelope
    for (itype i = 0; i < min(query->size(), subject->size()); ++i) {
        
        ftype cost = 0;
        
        if (onquery)
            cost = L->at(i) > subject->at(i) ? L->at(i)-subject->at(i) :
                   (U->at(i) < subject->at(i) ? subject->at(i)-U->at(i) : 0);
        else
            cost = L->at(i) > query->at(i) ? L->at(i)-query->at(i) :
                   (U->at(i) < query->at(i) ? query->at(i)-U->at(i) : 0);
    
        if (squared)
            cost *= cost;
        
        penalty += cost;
    }
    
    // free memory if needed
    if (!precalculated) {
        delete L;
        delete U;
    } else {
        if (onquery)
            delete query;
        else
            delete subject;
    }
    
    return penalty;
}

///////////////////////////////////////////////////////////////////////////////
// wrappers for lower bound - related methods
///////////////////////////////////////////////////////////////////////////////


int lb_envelope (std::vector<float> *series, unsigned int w, 
              std::vector<float> *L, std::vector<float> *U) {

    return meta_envelope<unsigned int, float> (series, w, L, U);
}

float lb_keogh_onQuery(std::vector<float> *query, std::vector<float> *subject,
                       unsigned int w, bool squared) {
    
    if (squared)
        return meta_lb_keogh<unsigned int, float, false, true, true> 
              (query, subject, w, 0, 0);
    else
        return meta_lb_keogh<unsigned int, float, false, true, false> 
              (query, subject, w,0 , 0);
}

float lb_keogh_onSubject(std::vector<float> *query, std::vector<float> *subject,
                         unsigned int w, bool squared) {
    
    if (squared)
        return meta_lb_keogh<unsigned int, float, false, false, true> 
              (query, subject, w, 0, 0);
    else
        return meta_lb_keogh<unsigned int, float, false, false, false> 
              (query, subject, w, 0, 0);
}

float lb_keogh_onEnvelope(std::vector<float> *series, std::vector<float> *L,
                          std::vector<float> *U, unsigned int w, bool squared) {
    if(squared)
        return meta_lb_keogh<unsigned int, float, true, true, true>
            (0, series, w, L, U);
    else
        return meta_lb_keogh<unsigned int, float, true, true, false>
            (0, series, w, L, U);
}


