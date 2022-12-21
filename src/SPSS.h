//
// Created by enrico on 21/12/22.
//

#ifndef USTAR_SPSS_H
#define USTAR_SPSS_H

#include "DBG.h"

using namespace std;

class SPSS{
    DBG *dbg;
    const vector<node_t> *nodes;

public:
    explicit SPSS(DBG *dbg);

    // void test();

    void simpler_test();
};

#endif //USTAR_SPSS_H
