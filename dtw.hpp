#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

float dist_euclidean(std::vector<float> *N, std::vector<float> *H);
float dist_manhatten(std::vector<float> *N, std::vector<float> *H);
float dist_dtw(std::vector<float> *N, std::vector<float> *H, bool squared);
float dist_dtw_backtrace(std::vector<float> *N, std::vector<float> *H, std::vector<std::pair<unsigned int, unsigned int> >* path, bool squared);
float dist_cdtw(std::vector<float> *N, std::vector<float> *H, unsigned int w ,bool squared);
float dist_cdtw_backtrace(std::vector<float> *N, std::vector<float> *H, unsigned int w, std::vector<std::pair<unsigned int, unsigned int> >* path, bool squared);

