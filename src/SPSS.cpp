//
// Created by enrico on 21/12/22.
//

#include <iostream>
#include <stack>
#include <fstream>
#include <algorithm>
#include "SPSS.h"

SPSS::SPSS(DBG *dbg, Sorter *sorter, bool duplicates ,bool debug){
    this->debug = debug;
    this->duplicates = duplicates;

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

void SPSS::extends(vector<node_idx_t> &path_nodes, vector<bool> &path_forwards) {
    vector<bool> dummy;
    vector<node_idx_t> to_nodes_e;
    vector<bool> to_forwards_e;

    vector<node_idx_t> to_nodes; vector<bool> to_forwards; vector<bool> forwards;

    // add the seed
    node_idx_t seed = path_nodes.back();
    bool forward = path_forwards.back();
    visited.at(seed) = true;

    // extension vars
    bool to_forward;
    size_t node, successor;

    // ----- extending -----
    node = seed;
    while (true) {
        // explore
        dbg->get_consistent_nodes_from(node, forward, to_nodes, to_forwards, visited);

        if (to_nodes.empty()) {
            if(!duplicates) break;

            // check visited nodes
            dbg->get_consistent_nodes_from(node, forward, to_nodes, to_forwards, dummy);
            if (to_nodes.empty()) break; // dead end

            // is there a node with non-visited neighbours?
            bool found = false;  bool dummy_f = false;
            int best = INT32_MAX; int conn = 0;
            node_idx_t best_node = to_nodes[0]; bool best_forward = to_forwards[0];
            vector<node_idx_t> best_to_nodes_e; vector<bool> best_to_forwards_e;
            for(int i = 0; i < to_nodes.size(); i++){
                node = to_nodes[i];
                forward = to_forwards[i];

                dbg->get_consistent_nodes_from(node, forward, to_nodes_e, to_forwards_e, visited);
                if(!to_nodes_e.empty()){ // node found!
                    found = true;
                    sorter->next_successor(node, forward, to_nodes_e, to_forwards_e, dummy_f, conn);
                    if(conn < best){
                        best = conn;
                        best_node = node;
                        best_forward = forward;
                        best_to_nodes_e = to_nodes_e;
                        best_to_forwards_e = to_forwards_e;
                        if(best == 1) break;
                    }
                }
            }
            if(!found)
                break; // no non-visited node
            else{
                // current node and neighbours
                node = best_node;
                forward = best_forward;
                to_nodes = best_to_nodes_e;
                to_forwards = best_to_forwards_e;
                // append to the path
                path_nodes.push_back(node);
                path_forwards.push_back(forward);
            }
        }
        successor = sorter->next_successor(node, forward, to_nodes, to_forwards, to_forward);

        // update
        node = successor;
        forward = to_forward;
        visited.at(node) = true;
        path_nodes.push_back(node);
        path_forwards.push_back(forward);
    }

    if(debug && !dbg->check_path_consistency(path_nodes, path_forwards))
        cerr << "Inconsistent path!\n";
}

void reverse_path(vector<node_idx_t> &path_nodes, vector<bool> &path_forwards){
    std::reverse(path_nodes.begin(), path_nodes.end());
    std::reverse(path_forwards.begin(), path_forwards.end());
    for(auto && forward : path_forwards)
        forward = (!forward);
}

void SPSS::compute_path_cover(bool two_way) {
    // reset path cover if already computed
    path_cover_nodes.clear();
    path_cover_forwards.clear();
    std::fill(visited.begin(), visited.end(), false);

    vector<node_idx_t> path_nodes; vector<bool> path_forwards;
    while(sorter->has_seed()){
        path_nodes.clear(); path_forwards.clear();

        path_nodes.push_back(sorter->next_seed());
        path_forwards.push_back(true);

        // forward extending
        extends(path_nodes, path_forwards);
        if(two_way) {
            if(debug) cout << "Doing two-way...\n";
            // backward extending
            reverse_path(path_nodes, path_forwards);
            extends(path_nodes, path_forwards);
            if(debug) cout << "Done.\n";
        }

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

    cout << "\n";
    cout << "Simplitigs stats:\n";
    cout << "   number of simplitigs (NS):              " << n_simplitigs << "\n";
    cout << "   cumulative length (CL):                 " << c_length << "\n";
    cout << "   nucleotide per kmer:                    " << (double) c_length / (double) dbg->get_n_kmers() << "\n";
    cout << "   average simplitigs length:              " << (double) c_length / (double) n_simplitigs << "\n";
    cout << "\n";
}

size_t SPSS::get_score() {
    return path_cover_nodes.size();
}
