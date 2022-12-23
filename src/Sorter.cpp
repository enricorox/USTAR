//
// Created by enrico on 22/12/22.
//

#include <iostream>
#include "Sorter.h"

Sorter::Sorter() {
    seed_order.resize(nodes->size());
    for(size_t i = 0; i < seed_order.size(); i++)
        seed_order[i] = i;

    seed_index = 0;
}

void Sorter::init(const vector<node_t> *nodes, const vector<bool> *visited){
    this->visited = visited;
    this->nodes = nodes;
}

size_t Sorter::next_seed() {
    if(!has_seed()){
        cerr << "Must check if a seed exist first!" << endl;
        exit(EXIT_FAILURE);
    }
    return seed_order.at(seed_index);
}

bool Sorter::has_seed() {
    for(; seed_index < seed_order.size(); seed_index++)
        if(!visited->at(seed_index))
            break;
    return seed_index < seed_order.size();
}
