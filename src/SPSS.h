//
// Created by enrico on 21/12/22.
//

#ifndef USTAR_SPSS_H
#define USTAR_SPSS_H

#include "DBG.h"

using namespace std;

class SPSS{
    DBG *dbg;

public:
    explicit SPSS(DBG *dbg);

    void simpler_test(const vector<size_t> &nodes, const vector<bool> &forwards);
};

#endif //USTAR_SPSS_H
