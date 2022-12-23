//
// Created by enrico on 22/12/22.
//

#include <iostream>
#include "Sorter.h"

Sorter::Sorter() {

}

void Sorter::init(const vector<node_t> *dbg_nodes, const vector<bool> *spss_visited){
    this->visited = spss_visited;
    this->nodes = dbg_nodes;

    seed_order.resize(dbg_nodes->size());
    for(size_t i = 0; i < seed_order.size(); i++)
        seed_order[i] = i;

    seed_index = 0;
}

size_t Sorter::next_seed() {
    if(!has_seed()){
        cerr << "next_seed(): Must check if a seed exists first!" << endl;
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
