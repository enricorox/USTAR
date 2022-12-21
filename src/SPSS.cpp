//
// Created by enrico on 21/12/22.
//

#include <iostream>
#include "SPSS.h"

SPSS::SPSS(DBG *dbg){
    this->dbg = dbg;
    nodes = dbg->get_nodes();
}

/*
void SPSS::test() {
    vector<const node_t *> path_nodes;
    vector<const arcs_t *> path_arcs;

    const node_t *node = &nodes->at(7);
    path_nodes.push_back(node);
    path_arcs.push_back(&node->arcs.at(0));

    node = &nodes->at(node->arcs.at(0).successor);
    path_nodes.push_back(node);
    path_arcs.push_back(&node->arcs.at(2));
    node = &nodes->at(node->arcs.at(2).successor);
    path_nodes.push_back(node);

    string contig = dbg->spell(path_nodes, path_arcs);
    cout << contig << " (" << contig.length() << ")" << endl;
}
*/
void SPSS::simpler_test() {
    vector<size_t> nodes = {7, 1578, 101636};
    vector<bool> forwards = {true, true, true};

    string contig = dbg->spell(nodes, forwards);
    cout << contig << " (" << contig.length() << ")" << endl;

    vector<uint32_t> counts;
    dbg->get_counts(nodes, forwards, counts);
    for(auto count : counts)
        cout << count << " ";
    cout << endl;
}
