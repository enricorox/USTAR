//
// Created by enrico on 21/12/22.
//

#ifndef USTAR_SPSS_H
#define USTAR_SPSS_H

#include "DBG.h"

using namespace std;

class SPSS{
    DBG *dbg;
    uint32_t n_nodes;

    // simplitigs
    vector<vector<size_t>> simplitigs_path_nodes;
    vector<vector<bool>> simplitigs_path_forwards;

    // saturated nodes
    vector<bool> saturated;

public:
    explicit SPSS(DBG *dbg);

    void test();
};

#endif //USTAR_SPSS_H
