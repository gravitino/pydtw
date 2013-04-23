#include <cmath>
#include <vector>
#include <iostream>

float euclidean(std::vector<float> *N, std::vector<float> *H);
float manhatten(std::vector<float> *N, std::vector<float> *H);
float dtw(std::vector<float> *N, std::vector<float> *H, bool squared=false);
float cdtw(std::vector<float> *N, std::vector<float> *H, unsigned int w, bool squared=false);
