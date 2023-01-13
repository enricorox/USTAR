//
// Created by enrico on 24/12/22.
//

#include <valarray>
#include <random>
#include "commons.h"

size_t get_rand(const size_t& a) {
    static std::random_device randDev;
    static std::mt19937 twister(randDev());
    static std::uniform_int_distribution<size_t> dist;

    dist.param(std::uniform_int_distribution<size_t>::param_type(0, a));
    return dist(twister);
}


uint32_t d(uint32_t a, uint32_t b){
    return abs((int)a - (int)b);
}

double d(double a, double b){
    return abs(a - b);
}

uint32_t median(vector<uint32_t> v){
    size_t n = v.size() / 2;
    // find the middle element
    nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}