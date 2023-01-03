//
// Created by enrico on 24/12/22.
//

#ifndef USTAR_COMMONS_H
#define USTAR_COMMONS_H

#include <cstdint>
#include <map>
#include <vector>

using namespace std;

template<typename T>
string inv_map(const map<string, T> &m, const T &name){
    for(auto &p : m)
        if(p.second == name)
            return p.first;
    return "?";
}

size_t get_rand(const size_t& a);

uint32_t d(uint32_t a, uint32_t b);

uint32_t median(vector<uint32_t> v);

#endif //USTAR_COMMONS_H
