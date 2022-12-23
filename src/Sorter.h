//
// Created by enrico on 22/12/22.
//

#ifndef USTAR_SORTER_H
#define USTAR_SORTER_H

#include <vector>
#include "DBG.h"

using namespace std;

class Sorter{
    const vector<bool> *visited{};
    const vector<node_t> *nodes{};

    vector<size_t> seed_order;
    size_t seed_index;

public:
    Sorter();

    void init(const vector<node_t> *dbg_nodes, const vector<bool> *spss_visited);

    bool has_seed();

    size_t next_seed();

    size_t next_successor(size_t node, bool forward, vector<size_t> to_nodes);
};

#endif //USTAR_SORTER_H
