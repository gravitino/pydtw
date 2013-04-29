#include <cmath>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) < (y) ? (y) : (x))
#define min3(x, y, z) (min(x, min(y, z)))
#define argmin3(x, y, z) ((x < y && x < z) ? 0 : (y < z ? 1 : 2))
#define sgn(x) ((x) < 0 ? -1 : 1)
#define abs(x) ((x)*sgn(x))
#define unused(x) ((void) x)


// distance measures

float dist_euclidean(std::vector<float> *N, std::vector<float> *H);
float dist_manhatten(std::vector<float> *N, std::vector<float> *H);
float dist_dtw(std::vector<float> *N, std::vector<float> *H, bool squared);
float dist_dtw_backtrace(std::vector<float> *N, std::vector<float> *H, std::vector<std::pair<unsigned int, unsigned int> >* path, bool squared);
float dist_cdtw(std::vector<float> *N, std::vector<float> *H, unsigned int w ,bool squared);
float dist_cdtw_backtrace(std::vector<float> *N, std::vector<float> *H, unsigned int w, std::vector<std::pair<unsigned int, unsigned int> >* path, bool squared);

// lower bounding techniques

int lb_envelope (std::vector<float> *series, unsigned int w, std::vector<float> *L, std::vector<float> *U);
float lb_keogh_onQuery(std::vector<float> *query, std::vector<float> *subject, unsigned int w, bool squared);
float lb_keogh_onSubject(std::vector<float> *query, std::vector<float> *subject, unsigned int w, bool squared);
float lb_keogh_onEnvelope(std::vector<float> *series, std::vector<float> *L, std::vector<float> *U, unsigned int w, bool squared);
