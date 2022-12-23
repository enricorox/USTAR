//
// Created by enrico on 22/12/22.
//

#include <iostream>
#include "Sorter.h"

Sorter::Sorter(const vector<sorting_method_t> &sorting_methods) {
    this->sorting_methods = sorting_methods;

}

void Sorter::init(const vector<node_t> *dbg_nodes, const vector<bool> *spss_visited){
    this->visited = spss_visited;
    this->nodes = dbg_nodes;

    seed_order.resize(dbg_nodes->size());
    for(size_t i = 0; i < seed_order.size(); i++)
        seed_order[i] = i;

    seed_index = 0;

    // sort!
    for(auto &sorting_method : sorting_methods)
        switch(sorting_method){
        case DEFAULT:
            // do nothing
            break;
        default:
            cerr << "init(): unknown sorting method!" << endl;
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

size_t Sorter::seed_successor(size_t seed, vector<bool> forwards, vector<size_t> to_nodes, vector<bool> to_forwards,
                              bool &forward, bool &to_forward) {
    size_t successor;
    forward = forwards[0];
    to_forward = to_forwards[0];
    successor = to_nodes[0];
    // scan all the sorting methods
    for(auto &sorting_method : sorting_methods)
        switch(sorting_method){
            case DEFAULT: // choose always the first
                // do nothing, it's before the cycle
                break;
            default:
                // do nothing
                break;
        }
    return successor;
}

size_t Sorter::next_successor(size_t node, vector<size_t> to_nodes, vector<bool> to_forwards, bool &to_forward) {
    // scan all the sorting methods
    size_t successor;
    to_forward = to_forwards[0];
    successor = to_nodes[0];

    for(auto &sorting_method : sorting_methods)
        switch(sorting_method){
            case DEFAULT:
                // do nothing: it's before the cycle
                break;
            default:
                // do nothing
                break;
        }
    return successor;
}
