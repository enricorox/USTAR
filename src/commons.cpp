//
// Created by enrico on 24/12/22.
//

#include <valarray>
#include "commons.h"

template<typename T>
string inv_map(const map<string, T> &m, const T &name){
    for(auto &p : m)
        if(p.second == name)
            return p.first;
    return "?";
}

uint32_t d(uint32_t a, uint32_t b){
    return abs((int)a - (int)b);
}

uint32_t median(vector<uint32_t> v){
    size_t n = v.size() / 2;
    // find the middle element
    nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}