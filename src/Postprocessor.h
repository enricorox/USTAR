//
// Created by enrico on 31/12/22.
//

#ifndef USTAR_POSTPROCESSOR_H
#define USTAR_POSTPROCESSOR_H

#include "consts.h"
#include <vector>

using namespace std;

class Postprocessor{
    vector<size_t> order;
    vector<bool> flips;

    size_t idx;

public:
    Postprocessor(const vector<string> &simplitigs, const vector<vector<uint32_t>> &counts);

    bool has_next();

    uint32_t next();
};

#endif //USTAR_POSTPROCESSOR_H
