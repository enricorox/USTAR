//
// Created by enrico on 22/12/22.
//

#include <iostream>
#include <algorithm>
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
        case seeding_method_t::FIRST:
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

size_t Sorter::seed_successor(node_idx_t seed, vector<bool> &forwards, vector<node_idx_t> &to_nodes, vector<bool> &to_forwards,
                              bool &forward, bool &to_forward) {

    vector<node_idx_t> order; order.reserve(to_nodes.size());
    for(node_idx_t i = 0; i < to_nodes.size(); i++)
        order.push_back(i); // default order

    switch(extending_method){
        case extending_method_t::FIRST: // choose always the first
            // do nothing, it's before the cycle
            break;
        case extending_method_t::SIMILAR_ABUNDANCE: {
            auto lambda = [this](node_idx_t a, node_idx_t b) {
                return nodes->at(a).abundances.at(0) < nodes->at(b).abundances.at(0);
            };
            sort(order.begin(), order.end(), lambda);
            }
            break;
        default:
            cerr << "seed_successor(): unknown extending method!" << endl;
            exit(EXIT_FAILURE);
    }
    if(!forwards.empty()) // is there a "dummy" forward?
        forward = forwards[order[0]];
    to_forward = to_forwards[order[0]];
    size_t successor = to_nodes[order[0]];
    return successor;
}

size_t Sorter::next_successor(node_idx_t node, vector<node_idx_t> &to_nodes, vector<bool> &to_forwards, bool &to_forward) {
    vector<bool> dummy_forwards;
    bool dummy_forward;
    return seed_successor(node, dummy_forwards, to_nodes, to_forwards, dummy_forward, to_forward);
}
