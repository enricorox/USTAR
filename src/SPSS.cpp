//
// Created by enrico on 21/12/22.
//

#include <iostream>
#include <stack>
#include <fstream>
#include "SPSS.h"

SPSS::SPSS(DBG *dbg){
    this->dbg = dbg;
    n_nodes = dbg->get_n_nodes();
    saturated.resize(n_nodes);
    fill(saturated.begin(), saturated.end(), false);
}

void SPSS::extends(uint32_t seed, vector<size_t> &path_nodes, vector<bool> &path_forwards, bool two_way) {
    path_nodes.clear(); path_forwards.clear();

    vector<size_t> to_nodes; vector<bool> to_forwards; vector<bool> forwards;
    deque<size_t> path_nodes_d; deque<bool> path_forwards_d;

    // forward extending
    dbg->get_nodes_from(seed, forwards, to_nodes, to_forwards, saturated);
    uint32_t node = seed;

    path_nodes_d.push_back(node);
    saturated.at(node) = true;
    if(to_nodes.empty()){
        path_nodes.push_back(seed);
        path_forwards.push_back(true);
        return;
    }
    // choose arbitrarily the first arc (if any!)
    bool seed_forward = forwards.at(0);
    bool forward = seed_forward;
    path_forwards_d.push_back(forward);
    while(true){
        dbg->get_consistent_nodes_from(node, forward, to_nodes, to_forwards, saturated);
        if(to_nodes.empty())
            break;
        node = to_nodes.at(0);
        forward = to_forwards.at(0);
        saturated.at(node) = true;
        path_nodes_d.push_back(node);
        path_forwards_d.push_back(forward);
    }

    // backward extending
    if(two_way) {
        node = seed;
        forward = (!seed_forward);
        while (true) {
            dbg->get_consistent_nodes_from(node, forward, to_nodes, to_forwards, saturated);
            if (to_nodes.empty())
                break;
            node = to_nodes.at(0);
            forward = to_forwards.at(0);
            saturated.at(node) = true;
            path_nodes_d.push_front(node);
            path_forwards_d.push_front(!forward);
        }
    }

    for(size_t i = 0; i < path_forwards_d.size(); i++){
        path_nodes.push_back(path_nodes_d.at(i));
        path_forwards.push_back(path_forwards_d.at(i));
    }
    if(!dbg->check_path_consistency(path_nodes, path_forwards))
        cerr << "Inconsistent path!\n";
}

void SPSS::extract_simplitigs() {
    // reset simplitigs if already computed
    simplitigs_path_nodes.clear();
    simplitigs_path_forwards.clear();

    vector<size_t> path_nodes; vector<bool> path_forwards;
    for(size_t seed = 0; seed < n_nodes; seed++){
        if(saturated[seed]) continue;
        extends(seed, path_nodes, path_forwards, true);
        simplitigs_path_nodes.push_back(path_nodes);
        simplitigs_path_forwards.push_back(path_forwards);
    }
    n_simplitigs = simplitigs_path_nodes.size();
}

void SPSS::to_fasta_file(const string &file_name) {
    if(n_simplitigs == 0){
        cerr << "Need to extract simplitigs first!" << endl;
        exit(EXIT_FAILURE);
    }
    ofstream fasta;
    fasta.open(file_name);
    for(size_t i = 0; i < n_simplitigs; i++){
        fasta << ">\n";
        fasta << dbg->spell(simplitigs_path_nodes[i], simplitigs_path_forwards[i]) << "\n";
    }
    fasta.close();
}

void SPSS::to_counts_file(const string &file_name) {
    ofstream counts_file;
    counts_file.open(file_name);
    vector<uint32_t> counts;
    for(size_t i = 0; i < n_simplitigs; i++){
        counts.clear();
        dbg->get_counts(simplitigs_path_nodes[i], simplitigs_path_forwards[i], counts);
        for(auto c : counts)
            counts_file << c << "\n";
    }
    counts_file.close();
}

void SPSS::print_info(){
    if(n_simplitigs == 0){
        cerr << "Need to extract simplitigs first!" << endl;
        exit(EXIT_FAILURE);
    }
    size_t c_length = 0;
    for(auto &simplitig : simplitigs_path_nodes){
        c_length += dbg->get_node(simplitig.at(0)).length;
        for(size_t i = 1; i < simplitig.size(); i++)
            c_length += dbg->get_node(simplitig.at(i)).length - dbg->get_kmer_size() + 1;
    }

    cout << "Simplitigs info:" << endl;
    cout << "\tnumber of simplitigs: NS = " << n_simplitigs << endl;
    cout << "\tcumulative length: CL = " << c_length << endl;
    cout << "\taverage simplitigs length: " << (double) c_length / (double) n_simplitigs << endl;
    cout << endl;
}
