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
    visited.resize(n_nodes);
    fill(visited.begin(), visited.end(), false);

    this->sorter = sorter;
    sorter->init(dbg->get_nodes(), &visited);
}

void SPSS::extends(size_t seed, vector<size_t> &path_nodes, vector<bool> &path_forwards, bool two_way) {
    path_nodes.clear(); path_forwards.clear();

    vector<size_t> to_nodes; vector<bool> to_forwards; vector<bool> forwards;
    deque<size_t> path_nodes_d; deque<bool> path_forwards_d;

    // forward extending
    dbg->get_nodes_from(seed, forwards, to_nodes, to_forwards, visited);
    uint32_t node = seed;

    path_nodes_d.push_back(node);
    visited.at(node) = true;
    if(to_nodes.empty()){
        path_nodes.push_back(seed);
        path_forwards.push_back(true);
        return;
    }
    // choose arbitrarily the first arc (if any!)
    bool seed_forward = forwards[0];
    bool forward = seed_forward;
    path_forwards_d.push_back(forward);
    while(true){
        dbg->get_consistent_nodes_from(node, forward, to_nodes, to_forwards, visited);
        if(to_nodes.empty())
            break;
        node = to_nodes[0];
        forward = to_forwards[0];
        visited.at(node) = true;
        path_nodes_d.push_back(node);
        path_forwards_d.push_back(forward);
    }

    // backward extending
    if(two_way) {
        node = seed;
        forward = (!seed_forward);
        while (true) {
            dbg->get_consistent_nodes_from(node, forward, to_nodes, to_forwards, visited);
            if (to_nodes.empty())
                break;
            node = to_nodes[0];
            forward = to_forwards[0];
            visited.at(node) = true;
            path_nodes_d.push_front(node);
            path_forwards_d.push_front(!forward);
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

    vector<size_t> path_nodes; vector<bool> path_forwards;
    for(size_t seed = 0; seed < n_nodes; seed++){
        if(visited[seed]) continue;
        extends(seed, path_nodes, path_forwards, two_way);
        path_cover_nodes.push_back(path_nodes);
        path_cover_forwards.push_back(path_forwards);
    }
}

void SPSS::extract_simplitigs_and_counts(){
    n_simplitigs = path_cover_nodes.size();
    if(n_simplitigs == 0){
        cerr << "Need to compute a path cover first!" << endl;
        exit(EXIT_FAILURE);
    }

    vector<size_t> simplitig_counts;
    for(size_t i = 0; i < n_simplitigs; i++){
        // extract simplitigs
        simplitigs.push_back(dbg->spell(path_cover_nodes[i], path_cover_forwards[i]));

        // extract counts
        simplitig_counts.clear();
        dbg->get_counts(path_cover_nodes[i], path_cover_forwards[i], simplitig_counts);
        counts.push_back(simplitig_counts);
    }
}

void SPSS::to_fasta_file(const string &file_name) {
    if(n_simplitigs == 0){
        cerr << "Need to extract simplitigs first!" << endl;
        exit(EXIT_FAILURE);
    }
    ofstream fasta;
    fasta.open(file_name);
    for(auto &simplitig : simplitigs){
        fasta << ">\n";
        fasta << simplitig << "\n";
    }
    fasta.close();
}

void SPSS::to_counts_file(const string &file_name) {
    ofstream counts_file;
    counts_file.open(file_name);
    for(const auto &simplitig_counts : counts){
        for(auto c : simplitig_counts)
            counts_file << c << (debug?" ":"\n");
        if(debug) counts_file << "\n";
    }
    counts_file.close();
}

void SPSS::print_info(){
    if(n_simplitigs == 0){
        cerr << "Need to extract simplitigs first!" << endl;
        exit(EXIT_FAILURE);
    }
    size_t c_length = 0;
    for(auto &simplitig : path_cover_nodes){
        c_length += dbg->get_node(simplitig[0]).length;
        for(size_t i = 1; i < simplitig.size(); i++)
            c_length += dbg->get_node(simplitig[i]).length - (dbg->get_kmer_size() - 1);
    }

    cout << "Simplitigs info:\n";
    cout << "\tnumber of simplitigs: NS = " << n_simplitigs << "\n";
    cout << "\tcumulative length: CL = " << c_length << "\n";
    cout << "\taverage simplitigs length: " << (double) c_length / (double) n_simplitigs << "\n";
    cout << endl;
}
