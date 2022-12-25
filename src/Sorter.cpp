//
// Created by enrico on 22/12/22.
//

#include <iostream>
#include <algorithm>
#include "Sorter.h"
#include "algos.h"

Sorter::Sorter(seeding_method_t sorting_methods, extending_method_t extending_method, bool debug) {
    this->seeding_method = sorting_methods;
    this->extending_method = extending_method;
    this->debug = debug;

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
        case seeding_method_t::LOWER_MEDIAN_ABUNDANCE: {
            auto lambda = [this](size_t a, size_t b) {
                return nodes->at(a).median_abundance < nodes->at(b).median_abundance;
            };
            if(debug)
                cout << "Sorting the seeds..." << endl;
            sort(seed_order.begin(), seed_order.end(), lambda);
            if(debug)
                cout << "Done." << endl;
            }
            break;
        case seeding_method_t::FIRST:
            // no break here
        case seeding_method_t::SIMILAR_ABUNDANCE:
            // do nothing
            break;
        default:
            cerr << "init(): unknown seeding method!" << endl;
            exit(EXIT_FAILURE);
    }
}

size_t Sorter::next_seed() {
    if(!has_seed()){
        cerr << "next_seed(): No seed available!" << endl;
        exit(EXIT_FAILURE);
    }
    return seed_order[seed_index];
}

bool Sorter::has_seed() {
    for(; seed_index < seed_order.size(); seed_index++)
        if(!(*visited)[seed_order[seed_index]])
            break;
    return seed_index < seed_order.size();
}

size_t Sorter::seed_successor(node_idx_t seed, vector<bool> &forwards, vector<node_idx_t> &to_nodes, vector<bool> &to_forwards,
                              bool &forward, bool &to_forward) {

    size_t best = 0;
    switch(extending_method){
        case extending_method_t::FIRST: // choose always the first
            // do nothing, it's before the cycle
            break;
        case extending_method_t::SIMILAR_ABUNDANCE: {
                uint32_t best_value = 999999999;
                for(size_t i = 0; i < to_nodes.size(); i++){
                    uint32_t ab_seed = nodes->at(seed).abundances.back();
                    uint32_t ab_succ = nodes->at(to_nodes.at(i)).abundances.front();

                    if(!forwards.at(i))
                        ab_seed = nodes->at(seed).abundances.front();
                    if(!to_forwards.at(i))
                        ab_succ = nodes->at(to_nodes.at(i)).abundances.back();

                    // compute the distance
                    uint32_t diff = d(ab_seed, ab_succ);

                    if(diff == 0){ // same abundance!
                        best = i;
                        break;
                    }
                    if(diff < best_value){
                        best_value = diff;
                        best = i;
                    }
                }
            }
            break;
        default:
            cerr << "seed_successor(): unknown extending method!" << endl;
            exit(EXIT_FAILURE);
    }

    forward = forwards[best];
    to_forward = to_forwards[best];
    size_t successor = to_nodes[best];
    return successor;
}

size_t Sorter::next_successor(node_idx_t node, bool forward, vector<node_idx_t> &to_nodes, vector<bool> &to_forwards, bool &to_forward) {
    vector<bool> dummy_forwards;
    dummy_forwards.resize(to_nodes.size(), forward);
    bool dummy_forward;
    return seed_successor(node, dummy_forwards, to_nodes, to_forwards, dummy_forward, to_forward);
}
