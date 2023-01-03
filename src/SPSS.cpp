//
// Created by enrico on 21/12/22.
//

#include <iostream>
#include <stack>
#include <fstream>
#include "SPSS.h"

SPSS::SPSS(DBG *dbg, Sorter *sorter, bool debug){
    this->debug = debug;

    this->dbg = dbg;
    n_nodes = dbg->get_n_nodes();
    visited.resize(n_nodes, false);

    this->sorter = sorter;
    sorter->init(dbg->get_nodes(), &visited);
}

const vector<string> * SPSS::get_simplitigs(){
    if(simplitigs.empty()){
        cerr << "There are no simplitigs!" << endl;
        exit(EXIT_FAILURE);
    }
    return &simplitigs;
}

const vector<vector<uint32_t>> * SPSS::get_counts(){
    if(simplitigs.empty()){
        cerr << "There are no simplitigs!" << endl;
        exit(EXIT_FAILURE);
    }
    return &counts;
}

void SPSS::extends(node_idx_t seed, vector<node_idx_t> &path_nodes, vector<bool> &path_forwards, bool two_way) {
    path_nodes.clear(); path_forwards.clear();

    vector<node_idx_t> to_nodes; vector<bool> to_forwards; vector<bool> forwards;
    deque<node_idx_t> path_nodes_d; deque<bool> path_forwards_d;

    // add the seed
    path_nodes_d.push_back(seed);
    visited.at(seed) = true;

    // ----- forward extending -----
    // get all the unvisited nodes reachable from seed
    dbg->get_nodes_from(seed, forwards, to_nodes, to_forwards, visited);

    if(to_nodes.empty()){
        path_nodes.push_back(seed);
        path_forwards.push_back(true);
        return;
    }

    // set the orientation of the seed according to the sorter
    bool seed_forward;
    bool forward, to_forward;
    size_t node;
    size_t successor = sorter->seed_successor(seed, forwards, to_nodes, to_forwards, seed_forward, to_forward);
    path_forwards_d.push_back(seed_forward);

    while(true){
        // update
        node = successor;
        visited.at(node) = true;
        forward = to_forward;
        path_nodes_d.push_back(node);
        path_forwards_d.push_back(forward);

        // explore
        dbg->get_consistent_nodes_from(node, forward, to_nodes, to_forwards, visited);
        if(to_nodes.empty())
            break;
        successor = sorter->next_successor(node, forward, to_nodes, to_forwards, to_forward);
    }

    // ----- backward extending -----
    if(two_way) {
        node = seed;
        forward = (!seed_forward); // the other side
        while (true) {
            // explore
            dbg->get_consistent_nodes_from(node, forward, to_nodes, to_forwards, visited);
            if (to_nodes.empty())
                break;
            successor = sorter->next_successor(node, forward, to_nodes, to_forwards, to_forward);

            // update
            node = successor;
            forward = to_forward;
            visited.at(node) = true;
            path_nodes_d.push_front(node);
            path_forwards_d.push_front(!forward); // need to be visited backward!!!
        }
    }

    // convert to vector
    for(size_t i = 0; i < path_forwards_d.size(); i++){
        path_nodes.push_back(path_nodes_d[i]);
        path_forwards.push_back(path_forwards_d[i]);
    }

    if(debug && !dbg->check_path_consistency(path_nodes, path_forwards))
        cerr << "Inconsistent path!\n";
}

void SPSS::compute_path_cover(bool two_way) {
    // reset path cover if already computed
    path_cover_nodes.clear();
    path_cover_forwards.clear();
    std::fill(visited.begin(), visited.end(), false);

    vector<node_idx_t> path_nodes; vector<bool> path_forwards;
    while(sorter->has_seed()){
        extends(sorter->next_seed(), path_nodes, path_forwards, two_way);
        path_cover_nodes.push_back(path_nodes);
        path_cover_forwards.push_back(path_forwards);
    }
}

void SPSS::extract_simplitigs_and_counts(){
    n_simplitigs = path_cover_nodes.size();
    if(n_simplitigs == 0){
        cerr << "extract_simplitigs_and_counts(): Need to compute a path cover first!" << endl;
        exit(EXIT_FAILURE);
    }

    vector<uint32_t> simplitig_counts;
    counts.reserve(n_simplitigs);
    for(size_t i = 0; i < n_simplitigs; i++){
        // extract simplitigs
        simplitigs.push_back(dbg->spell(path_cover_nodes[i], path_cover_forwards[i]));

        // extract simplitigs_counts
        simplitig_counts.clear();
        dbg->get_counts(path_cover_nodes[i], path_cover_forwards[i], simplitig_counts);
        counts.push_back(simplitig_counts);
    }
}

void SPSS::print_stats(){
    if(n_simplitigs == 0){
        cerr << "print_stats(): Need to extract simplitigs first!" << endl;
        exit(EXIT_FAILURE);
    }
    size_t c_length = 0;
    for(auto &path : path_cover_nodes){
        c_length += dbg->get_node(path[0]).length;
        for(size_t i = 1; i < path.size(); i++)
            c_length += dbg->get_node(path[i]).length - (dbg->get_kmer_size() - 1);
    }

    cout << "\nSimplitigs stats:\n";
    cout << "   number of simplitigs (NS):              " << n_simplitigs << "\n";
    cout << "   cumulative length (CL):                 " << c_length << "\n";
    cout << "   nucleotide per kmer:                    " << (double) c_length / (double) dbg->get_n_kmers() << "\n";
    cout << "   average simplitigs length:              " << (double) c_length / (double) n_simplitigs << "\n";
    cout << "   average number of kmers per simplitig:  " << (double) dbg->get_n_kmers() / (double) n_simplitigs << "\n";
    cout << endl;
}

size_t SPSS::get_score() {
    return path_cover_nodes.size();
}
