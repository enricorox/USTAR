//
// Created by enrico on 22/12/22.
//

#ifndef USTAR_SORTER_H
#define USTAR_SORTER_H

#include <vector>
#include "DBG.h"

using namespace std;

enum sorting_method_t{
    DEFAULT
};

class Sorter{
    vector<sorting_method_t> sorting_methods;

    const vector<bool> *visited{};
    const vector<node_t> *nodes{};

    vector<size_t> seed_order;
    size_t seed_index = 0;

public:
    explicit Sorter(const vector<sorting_method_t> &sorting_methods={DEFAULT});

    void init(const vector<node_t> *dbg_nodes, const vector<bool> *spss_visited);

    bool has_seed();

    size_t next_seed();

    size_t seed_successor(size_t seed, vector<bool> forwards, vector<size_t> to_nodes, vector<bool> to_forwards, bool &forward, bool &to_forward);

    size_t next_successor(size_t node, vector<size_t> to_nodes, vector<bool> to_forwards, bool &to_forward);
};

#endif //USTAR_SORTER_H
