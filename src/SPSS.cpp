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

void SPSS::test() {
    vector<vector<size_t>> nodes_v = {{7, 157883, 101636}, {7, 157883, 82391}, {7, 157883, 114132}};
    vector<vector<bool>> forwards_v = {{true, true, true}, {true, true, true}, {true, true, true}};

    for(size_t i = 0; i < nodes_v.size(); i++) {
        auto &nodes = nodes_v.at(i);
        auto &forwards = forwards_v.at(i);

        if (dbg->check_path_consistency(nodes, forwards))
            cout << "Path is consistent!" << endl;
        else
            cout << "Path is NOT consistent!" << endl;

        string contig = dbg->spell(nodes, forwards);
        cout << contig << " (" << contig.length() << ")" << endl;

        vector<uint32_t> counts;
        dbg->get_counts(nodes, forwards, counts);
        for (auto count: counts)
            cout << count << " ";
        cout << endl;
    }
    cout << "extending node 18..." << endl;
    vector<size_t> path_nodes; vector<bool> path_forwards;
    extends(18, path_nodes, path_forwards);
    string contig = dbg->spell(path_nodes, path_forwards);
    cout << contig << " (" << contig.length() << ")" << endl;
}

void SPSS::extends(uint32_t seed, vector<size_t> &path_nodes, vector<bool> &path_forwards) {
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
    node = seed;
    forward = seed_forward;
    while(true){
        dbg->get_consistent_nodes_from(node, !forward, to_nodes, to_forwards, saturated);
        if(to_nodes.empty())
            break;
        node = to_nodes.at(0);
        forward = to_forwards.at(0);
        saturated.at(node) = true;
        path_nodes_d.push_front(node);
        path_forwards_d.push_front(forward);
    }

    for(size_t i = 0; i < path_forwards_d.size(); i++){
        path_nodes.push_back(path_nodes_d.at(i));
        path_forwards.push_back(path_forwards_d.at(i));
    }
}

void SPSS::extract_simplitigs() {
    vector<size_t> path_nodes; vector<bool> path_forwards;
    for(size_t seed = 0; seed < n_nodes; seed++){
        if(saturated[seed]) continue;
        extends(seed, path_nodes, path_forwards);
        simplitigs_path_nodes.push_back(path_nodes);
        simplitigs_path_forwards.push_back(path_forwards);
    }
    n_simplitigs = simplitigs_path_nodes.size();
}

void SPSS::to_fasta(const string &file_name) {
    if(n_simplitigs == 0){
        cerr << "Need to extract simplitigs first!" << endl;
        exit(EXIT_FAILURE);
    }
    ofstream fasta;
    fasta.open(file_name);
    for(size_t i = 0; i < n_simplitigs; i++){
        fasta << ">\n";
        fasta << dbg->spell(simplitigs_path_nodes.at(i), simplitigs_path_forwards.at(i)) << "\n";
    }
    fasta.close();
}