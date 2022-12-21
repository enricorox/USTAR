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

        // ------ parse edges ------
        // token = "L:-:23:+ L:-:104831:+  L:+:22:-"
        while(token != nullptr){
            edge_t edge{};
            char c1, c2;
            sscanf(token, "%*2c %c %*c %d %*c %c", &c1, &edge.successor, &c2); // L:-:23:+
            edge.forward = (c1 == '+');
            edge.to_forward = (c2 == '+');
            node.edges.push_back(edge);
            // next edges
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
    for(const auto &node : nodes) {
        n_edges += node.edges.size();
        n_kmers += node.abundances.size();
        sum_unitig_length += node.length;
    }
    avg_unitig_len = (double) sum_unitig_length / (double) nodes.size();
}

DBG::~DBG() = default;

void DBG::print_info() {
    cout << "Info for " << bcalm_file_name << ":" << endl;
    cout << "\tnumber of kmers: " << n_kmers << endl;
    cout << "\tnumber of nodes: " << nodes.size() << endl;
    cout << "\tnumber of edges: " << n_edges << endl;
    cout << "\taverage number of edges: " << (double) n_edges / (double) nodes.size() << endl;
    cout << "\taverage unitig length: " << avg_unitig_len << endl;
}

bool DBG::verify_overlaps() {
    for(const auto &node : nodes){
        for(const auto &edge : node.edges)
            if(!overlaps(node, edge))
                return false;
    }
    return true;
}

bool DBG::overlaps(const node_t &node, const edge_t &edge){
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
        for(auto &edge : node.edges)
            file << "L:" << (edge.forward?"+":"-") << ":" << edge.successor << ":" << (edge.to_forward?"+":"-") << " ";
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

string DBG::spell(const vector<node_t *> &path_nodes, const vector<edge_t *> &path_edges) {
    if(!check_path_consistency(path_nodes, path_edges))
        return "";

    string simplitig;
    if(path_edges.at(0)->forward)
        simplitig = path_nodes.at(0)->unitig;
    else
        simplitig = reverse_complement(path_nodes.at(0)->unitig);

    for(int i = 1; i < path_edges.size(); i++){

    }
    return simplitig;
}

bool DBG::check_path_consistency(const vector<node_t *> &path_nodes, const vector<edge_t *> &path_edges) {
    if(path_edges.size() != path_nodes.size() - 1)
        return false;

    bool last_forward = path_edges.at(0)->to_forward;
    for(int i = 1; i < path_edges.size(); i++){
        // check node-->edge-->node consistency
        if(&nodes.at(path_edges.at(i - 1)->successor) != path_nodes.at(i))
            return false;

        // check edge signs consistency
        if(last_forward != path_edges.at(i)->forward)
            return false;
        last_forward = path_edges.at(i)->to_forward;
    }
    return true;
}
