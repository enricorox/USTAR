//
// Created by enrico on 20/12/22.
//

#include <fstream>
#include <iostream>
#include <cstring>
#include "DBG.h"

void DBG::parse_bcalm_file() {
    ifstream bcalm_file;
    bcalm_file.open(bcalm_file_name);

    if(!bcalm_file.good()){
        cerr << "Can't access file " << bcalm_file_name << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    int serial;
    char dyn_line[MAX_LINE_LEN];
    while(getline(bcalm_file, line)){
        // make a new node
        node_t node;

        // check if line fits in buffer
        if(line.size() > MAX_LINE_LEN){
            cerr << "Lines must be smaller than " << MAX_LINE_LEN << " characters!" << endl;
            exit(EXIT_FAILURE);
        }

        // ------ parse line ------
        // >25 LN:i:32 ab:Z:14 12   L:-:23:+ L:-:104831:+  L:+:22:-

        // check consistency
        if(line[0] != '>'){
            cerr << "Bad formatted input file: no def-line found!" << endl;
            exit(EXIT_FAILURE);
        }

        // find BCALM2 sequence serial and unitig length
        // format: (ignore 1 char) (read 1 integer) (ignore 5 char) (read 1 integer) (read 1 string)
        sscanf(line.c_str(), "%*c %d %*5c %d %[^\n]s", &serial, &node.length, dyn_line);

        // check consistency
        if(serial != nodes.size()){
            cerr << "Bad formatted input file: lines must have progressive IDs!" << endl;
            exit(EXIT_FAILURE);
        }

        // ------ parse abundances ------
        // dyn_line = "ab:Z:14 12   L:-:23:+ L:-:104831:+  L:+:22:-"
        int sum = 0;
        char *token = strtok(dyn_line + 5, " ");
        do{
            int abundance = atoi(token);
            sum += abundance;
            node.abundances.push_back(abundance);
            token = strtok(nullptr, " ");
        }while(token != nullptr && token[0] != 'L');
        node.avg_abundance = sum / (double) node.abundances.size();

        // ------ parse arcs ------
        // token = "L:-:23:+ L:-:104831:+  L:+:22:-"
        while(token != nullptr){
            arcs_t edge{};
            char c1, c2;
            sscanf(token, "%*2c %c %*c %d %*c %c", &c1, &edge.successor, &c2); // L:-:23:+
            edge.forward = (c1 == '+');
            edge.to_forward = (c2 == '+');
            node.arcs.push_back(edge);
            // next arcs
            token = strtok(nullptr, " ");
        }

        // ------ parse sequence line ------
        // TTGAAGGTAACGGATGTTCTAGTTTTTTCTCTTT}
        getline(bcalm_file, line);

        // get the sequence
        node.unitig = line;

        // check consistency
        if((node.unitig.size() - kmer_size + 1) != node.abundances.size()){
            cerr << "Bad formatted input file: wrong number of abundances!" << endl;
            cerr << "Also make sure that kmer_size=" << kmer_size << endl;
            exit(EXIT_FAILURE);
        }

        // save the node
        nodes.push_back(node);
    }

    bcalm_file.close();
}

DBG::DBG(const string &bcalm_file_name, int kmer_size){
    this->bcalm_file_name = bcalm_file_name;
    this->kmer_size = kmer_size;

    // build the graph
    parse_bcalm_file();

    // compute graph parameters
    size_t sum_unitig_length = 0;
    double sum_abundances = 0;
    for(const auto &node : nodes) {
        n_arcs += node.arcs.size();
        n_kmers += node.abundances.size();
        sum_unitig_length += node.length;
        sum_abundances += node.avg_abundance * node.abundances.size();
    }
    avg_unitig_len = (double) sum_unitig_length / (double) nodes.size();
    avg_abundances = sum_abundances / n_kmers;
}

DBG::~DBG() = default;

void DBG::print_info() {
    cout << "Info for " << bcalm_file_name << ":" << endl;
    cout << "\tnumber of kmers: " << n_kmers << endl;
    cout << "\tnumber of nodes: " << nodes.size() << endl;
    cout << "\tnumber of arcs: " << n_arcs << endl;
    cout << "\taverage number of arcs: " << (double) n_arcs / (double) nodes.size() << endl;
    cout << "\taverage unitig length: " << avg_unitig_len << endl;
    cout << "\taverage abundances: " << avg_abundances << endl;
}

bool DBG::verify_overlaps() {
    for(const auto &node : nodes){
        for(const auto &edge : node.arcs)
            if(!overlaps(node, edge))
                return false;
    }
    return true;
}

bool DBG::overlaps(const node_t &node, const arcs_t &edge){
    string u1, u2;
    if (edge.forward) // + --> +/-
        u1 = node.unitig.substr(node.unitig.length() - kmer_size + 1);
    else // - --> +/-
        u1 = reverse_complement(node.unitig.substr(0, kmer_size - 1));

    if(edge.to_forward) // +/- --> +
        u2 = nodes[edge.successor].unitig.substr(0, kmer_size - 1);
    else  // +/- --> -
        u2 = reverse_complement(nodes[edge.successor].unitig.substr(nodes[edge.successor].unitig.length() - kmer_size + 1));

    return u1 == u2;
}

string DBG::reverse_complement(const string &s) {
    string rc(s);

    char c;
    for(int i = 0; i < s.length(); i++) {
        switch (s[i]) {
            case 'A':
            case 'a':
                c = 'T';
                break;
            case 'C':
            case 'c':
                c = 'G';
                break;
            case 'T':
            case 't':
                c = 'A';
                break;
            case 'G':
            case 'g':
                c = 'C';
                break;
        }
        rc[s.length() - 1 - i] = c;
    }
    return rc;
}

void DBG::to_bcalm_file(const string &file_name) {
    ofstream file;
    file.open(file_name);

    int id = 0;
    for(const auto &node : nodes){
        // >3 LN:i:33 ab:Z:2 2 3    L:+:138996:+
        // CAAAACCAGACATAATAAAAATACTAATTAATG
        file << ">" << id++ << " LN:i:" << node.length << " ab:Z:";
        for(auto &ab : node.abundances)
            file << ab << " ";
        for(auto &arcs : node.arcs)
            file << "L:" << (arcs.forward ? "+" : "-") << ":" << arcs.successor << ":" << (arcs.to_forward ? "+" : "-") << " ";
        file << "\n" << node.unitig << "\n";
    }

    file.close();
}

bool DBG::validate(){
    string test = "USTAR-test.fasta";
    to_bcalm_file(test);

    ifstream s1, s2;
    s1.open(bcalm_file_name);
    s2.open(test);

    string tok1, tok2;
    while(s1 >> tok1 && s2 >> tok2)
        if(tok1 != tok2) {
            cout << tok1 << " (" << tok1.length() <<") != " << tok2 << " (" << tok2.length() << ")" << endl;
            return false;
        }
    return true;
}

/*
string DBG::spell(const vector<const node_t *> &path_nodes, const vector<const arcs_t *> &path_arcs) {
    if(!check_path_consistency(path_nodes, path_arcs)) {
        cerr << "This path is not consistent!" << endl;
        exit(EXIT_FAILURE);
    }else
        cout << "Path is consistent" << endl;

    string simplitig;
    if(path_arcs.empty() || path_arcs.at(0)->forward)
        simplitig = path_nodes.at(0)->unitig;
    else
        simplitig = reverse_complement(path_nodes.at(0)->unitig);

    for(int i = 0; i < path_arcs.size(); i++){
        if(path_arcs.at(i)->to_forward)
            simplitig += path_nodes.at(i + 1)->unitig.substr(kmer_size - 1);
        else
            simplitig += reverse_complement(path_nodes.at(i + 1)->unitig.substr(kmer_size - 1));
    }
    return simplitig;
}
*/
bool DBG::check_path_consistency(const vector<const node_t *> &path_nodes, const vector<const arcs_t *> &path_arcs) {
    if(path_arcs.size() != path_nodes.size() - 1) {
        return false;
    }
    if(path_arcs.empty())
        return true;

    cout << "checking consistency..." << endl;
    bool last_forward = path_arcs.at(0)->to_forward;
    for(int i = 1; i < path_arcs.size(); i++){
        // check node-->edge-->node consistency
        if(&nodes.at(path_arcs.at(i - 1)->successor) != path_nodes.at(i))
            return false;

        // check edge signs consistency
        if(last_forward != path_arcs.at(i)->forward)
            return false;
        last_forward = path_arcs.at(i)->to_forward;
    }
    return true;
}

const vector<node_t> * DBG::get_nodes(){
    return &nodes;
}

void DBG::get_nodes_from(int node, vector<size_t> &to_nodes, vector<bool> &forwards) {
    to_nodes.clear();
    forwards.clear();
    for(auto &arc : nodes.at(node).arcs){
        to_nodes.push_back(arc.successor);
        forwards.push_back(arc.to_forward);
    }
}

void DBG::get_consistent_nodes_from(int node, bool forward, vector<size_t> &to_nodes, vector<bool> &forwards) {
    to_nodes.clear();
    forwards.clear();
    for(auto &arc : nodes.at(node).arcs){
        if(arc.forward == forward) {
            to_nodes.push_back(arc.successor);
            forwards.push_back(arc.to_forward);
        }
    }
}

string DBG::spell(const vector<size_t> &path_nodes, const vector<bool> &forwards) {
    if(path_nodes.size() != forwards.size()){
        cerr << "Inconsistent path!" << endl;
        exit(EXIT_FAILURE);
    }
    if(path_nodes.empty())
        return "";

    string simplitig;
    if(forwards.at(0))
        simplitig = nodes.at(path_nodes.at(0)).unitig;
    else
        simplitig = reverse_complement(nodes.at(path_nodes.at(0)).unitig);

    for(int i = 1; i < path_nodes.size(); i++){
        if(forwards.at(i))
            simplitig += nodes.at(path_nodes.at(i)).unitig.substr(kmer_size - 1);
        else
            simplitig += reverse_complement(nodes.at(path_nodes.at(i)).unitig.substr(kmer_size - 1));
    }

    return simplitig;
}

void DBG::get_counts(const vector<size_t> &path_nodes, const vector<bool> &forwards, vector<uint32_t> &counts) {
    for (int i = 0; i < path_nodes.size(); i++)
        if (forwards.at(i))
            for (size_t k = 0; k < nodes.at(path_nodes.at(i)).abundances.size(); k++)
                counts.push_back(nodes.at(path_nodes.at(i)).abundances.at(k));
        else
            for (size_t k = nodes.at(path_nodes.at(i)).abundances.size() - 1; k > -1; k--)
                counts.push_back(nodes.at(path_nodes.at(i)).abundances.at(k));
}


