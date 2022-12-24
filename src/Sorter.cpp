//
// Created by enrico on 22/12/22.
//

#include <iostream>
#include "Sorter.h"

Sorter::Sorter(seeding_method_t sorting_methods, extending_method_t extending_method) {
    this->seeding_method = sorting_methods;
    this->extending_method = extending_method;

}

void Sorter::init(const vector<node_t> *dbg_nodes, const vector<bool> *spss_visited){
    this->visited = spss_visited;
    this->nodes = dbg_nodes;

    seed_order.reserve(dbg_nodes->size());
    for(size_t i = 0; i < dbg_nodes->size(); i++)
        seed_order.push_back(i);

    seed_index = 0;

    // sort!
    switch(seeding_method){
        case seeding_method_t::DEFAULT_SEED:
        // do nothing
        break;
    default:
        cerr << "init(): unknown seeding method!" << endl;
        exit(EXIT_FAILURE);
    }
}

size_t Sorter::next_seed() {
    if(!has_seed()){
        cerr << "next_seed(): Must check if a seed exists first!" << endl;
        exit(EXIT_FAILURE);
    }
    return seed_order[seed_index];
}

bool Sorter::has_seed() {
    for(; seed_index < seed_order.size(); seed_index++)
        if(!(*visited)[seed_index])
            break;
    return seed_index < seed_order.size();
}

size_t Sorter::seed_successor(node_idx_t seed, vector<bool> forwards, vector<node_idx_t> to_nodes, vector<bool> to_forwards,
                              bool &forward, bool &to_forward) {
    size_t successor;
    forward = forwards[0];
    to_forward = to_forwards[0];
    successor = to_nodes[0];
    // scan all the sorting methods

    switch(extending_method){
        case extending_method_t::DEFAULT_EXTEND: // choose always the first
            // do nothing, it's before the cycle
            break;
        default:
            cerr << "seed_successor(): unknown extending method!" << endl;
            exit(EXIT_FAILURE);
    }
    return successor;
}

size_t Sorter::next_successor(node_idx_t node, vector<node_idx_t> to_nodes, vector<bool> to_forwards, bool &to_forward) {
    // scan all the sorting methods
    size_t successor;
    to_forward = to_forwards[0];
    successor = to_nodes[0];

    switch(extending_method){
        case extending_method_t::DEFAULT_EXTEND:
            // do nothing: it's before the cycle
            break;
        default:
            cerr << "next_successor(): unknown extending method!" << endl;
            exit(EXIT_FAILURE);
    }
    return successor;
}
