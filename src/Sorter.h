//
// Created by enrico on 22/12/22.
//

#ifndef USTAR_SORTER_H
#define USTAR_SORTER_H

#include <vector>
#include <random>
#include "DBG.h"
#include "consts.h"

using namespace std;

class Sorter{
    mt19937 random_generator;
    bool debug;
    seeding_method_t seeding_method;
    extending_method_t extending_method;

    const vector<bool> *visited{};
    const vector<node_t> *nodes{};

    vector<size_t> seed_order;
    size_t seed_index = 0;
    size_t last_node = 0;
    bool first_node = true;

public:
    Sorter(seeding_method_t seeding_method, extending_method_t extending_method, bool debug=false);

    void init(const vector<node_t> *dbg_nodes, const vector<bool> *spss_visited);

    bool has_seed();

    size_t next_seed();

    size_t next_successor(node_idx_t node, bool forward, vector<node_idx_t> &to_nodes, vector<bool> &to_forwards, bool &to_forward, int &best);
    size_t next_successor(node_idx_t node, bool forward, vector<node_idx_t> &to_nodes, vector<bool> &to_forwards, bool &to_forward);
};

#endif //USTAR_SORTER_H
