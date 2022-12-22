//
// Created by enrico on 20/12/22.
//

#include <fstream>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "DBG.h"

uint32_t static median(vector<uint32_t> v){
    size_t n = v.size() / 2;
    // find the middle element
    nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}

void DBG::parse_bcalm_file() {
    ifstream bcalm_file;
    bcalm_file.open(bcalm_file_name);

    if(!bcalm_file.good()){
        cerr << "Can't access file " << bcalm_file_name << endl;
        exit(EXIT_FAILURE);
    }

    // start parsing two line at a time
    string line;
    while(getline(bcalm_file, line)){
        // escape comments
        if(line[0] == '#')
            continue;

        size_t serial; // BCALM2 serial
        char dyn_line[MAX_LINE_LEN]; // line after id and length

        // make a new node
        node_t node;

        // check if line fits in dyn_line
        if(line.size() > MAX_LINE_LEN){
            cerr << "Lines must be smaller than " << MAX_LINE_LEN << " characters!" << endl;
            exit(EXIT_FAILURE);
        }

        // ------ parse line ------
        // >25 LN:i:32 ab:Z:14 12   L:-:23:+ L:-:104831:+  L:+:22:-

        // check consistency:
        // must have a def-line
        if(line[0] != '>'){
            cerr << "Bad formatted input file: no def-line found!" << endl;
            exit(EXIT_FAILURE);
        }

        // find BCALM2 sequence serial and unitig length
        // format: (ignore 1 char) (read 1 integer) (ignore 5 char) (read 1 integer) (read 1 string)
        sscanf(line.c_str(), "%*c %zd %*5c %d %[^\n]s", &serial, &node.length, dyn_line);

        // check consistency:
        // must have progressive IDs
        if(serial != nodes.size()){
            cerr << "Bad formatted input file: lines must have progressive IDs!" << endl;
            exit(EXIT_FAILURE);
        }

        // ------ parse abundances ------
        // dyn_line = "ab:Z:14 12   L:-:23:+ L:-:104831:+  L:+:22:-"
        uint32_t sum_abundance = 0;
        char *token = strtok(dyn_line + 5, " "); // tokenize abundances
        do{
            uint32_t abundance = atoi(token);
            sum_abundance += abundance;
            node.abundances.push_back(abundance);
            token = strtok(nullptr, " "); // next token
        }while(token != nullptr && token[0] != 'L');
        node.mean_abundance = sum_abundance / (double) node.abundances.size();
        node.median_abundance = median(node.abundances);

        // ------ parse arcs ------
        // token = "L:-:23:+ L:-:104831:+  L:+:22:-"
        while(token != nullptr){
            arc_t arc{};
            char s1, s2; // left and right signs
            sscanf(token, "%*2c %c %*c %d %*c %c", &s1, &arc.successor, &s2); // L:-:23:+
            arc.forward = (s1 == '+');
            arc.to_forward = (s2 == '+');
            node.arcs.push_back(arc);
            // next arcs
            token = strtok(nullptr, " ");
        }

        // ------ parse sequence line ------
        // TTGAAGGTAACGGATGTTCTAGTTTTTTCTCTTT}
        getline(bcalm_file, line);

        // get the sequence
        node.unitig = line;

        // check consistency:
        // there must be one count for each kmer!
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

DBG::DBG(const string &bcalm_file_name, uint32_t kmer_size){
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
        sum_abundances += node.mean_abundance * (double) node.abundances.size();
    }
    avg_unitig_len = (double) sum_unitig_length / (double) nodes.size();
    avg_abundances = sum_abundances / (double) n_kmers;
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
    cout << endl;
}

bool DBG::verify_overlaps() {
    for(const auto &node : nodes){
        for(const auto &arc : node.arcs)
            if(!overlaps(node, arc))
                return false;
    }
    return true;
}

bool DBG::overlaps(const node_t &node, const arc_t &arcs){
    string u1, u2;
    if (arcs.forward) // + --> +/-
        // last kmer_size - 1 characters
        u1 = node.unitig.substr(node.unitig.length() - kmer_size + 1);
    else // - --> +/-
        // first kmer_size - 1 characters reverse-complemented
        u1 = reverse_complement(node.unitig.substr(0, kmer_size - 1));

    if(arcs.to_forward) // +/- --> +
        // first kmer_size - 1 characters
        u2 = nodes[arcs.successor].unitig.substr(0, kmer_size - 1);
    else  // +/- --> -
        // last kmer_size - 1 characters reverse-complemented
        u2 = reverse_complement(nodes[arcs.successor].unitig.substr(nodes[arcs.successor].unitig.length() - kmer_size + 1));

    return u1 == u2;
}

string DBG::reverse_complement(const string &s) {
    string rc(s);

    for(size_t i = 0; i < s.length(); i++) {
        char c;
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
            default:
                cerr << "Unknown nucleotide!" << endl;
                exit(EXIT_FAILURE);
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
    string fasta_dbg = "unitigs.k"+ to_string(kmer_size) +".ustar.fa";
    to_bcalm_file(fasta_dbg);

    ifstream bcalm_dbg, this_dbg;
    bcalm_dbg.open(bcalm_file_name);
    this_dbg.open(fasta_dbg);

    string tok1, tok2;
    while(bcalm_dbg >> tok1 && this_dbg >> tok2)
        if(tok1 != tok2) {
            cerr << "Files differ here: " << tok1 << " != " <<  tok2.length() << endl;
            return false;
        }
    return true;
}

bool DBG::verify_input(){
    bool good = true;
    if (verify_overlaps())
        cout << "YES! DBG is an overlapping graph!" << endl;
    else {
        cout << "OOPS! DBG is NOT an overlapping graph" << endl;
        good = false;
    }
    if(validate())
        cout << "YES! DBG is the same as BCALM2 one!" << endl;
    else {
        cout << "OOPS! DBG is NOT the same as BCALM2 one!" << endl;
        good = false;
    }
    return good;
}

void DBG::get_nodes_from(uint32_t node, vector<bool> &forwards, vector<size_t> &to_nodes, vector<bool> &to_forwards, const vector<bool> &mask) {
    // vectors must be empty
    to_nodes.clear();
    to_forwards.clear();

    for(auto &arc : nodes.at(node).arcs){
        if(mask.at(arc.successor)) continue;
        to_nodes.push_back(arc.successor);
        forwards.push_back(arc.forward);
        to_forwards.push_back(arc.to_forward);
    }
}

void DBG::get_consistent_nodes_from(uint32_t node, bool forward, vector<size_t> &to_nodes, vector<bool> &to_forwards, const vector<bool> &mask) {
    // vectors must be empty
    to_nodes.clear();
    to_forwards.clear();

    for(auto &arc : nodes.at(node).arcs){
        if(mask.at(arc.successor)) continue;
        if(arc.forward == forward) { // consistent nodes only
            to_nodes.push_back(arc.successor);
            to_forwards.push_back(arc.to_forward);
        }
    }
}

string DBG::spell(const vector<size_t> &path_nodes, const vector<bool> &forwards) {
    if(path_nodes.size() != forwards.size()){
        cerr << "Inconsistent path!" << endl;
        exit(EXIT_FAILURE);
    }
    if(path_nodes.empty()) {
        cerr << "You're not allowed to spell an empty path!" << endl;
        exit(EXIT_FAILURE);
    }

    string contig;
    // first node as a seed
    if(forwards.at(0))
        contig = nodes.at(path_nodes.at(0)).unitig;
    else
        contig = reverse_complement(nodes.at(path_nodes.at(0)).unitig);

    // extend the seed
    for(size_t i = 1; i < path_nodes.size(); i++){
        if(forwards.at(i))
            contig += nodes.at(path_nodes.at(i)).unitig.substr(kmer_size - 1);
        else {
            string unitig = nodes.at(path_nodes.at(i)).unitig;
            size_t len = unitig.length() - kmer_size + 1;
            contig += reverse_complement(unitig.substr(0, len));
        }
    }

    return contig;
}

void DBG::get_counts(const vector<size_t> &path_nodes, const vector<bool> &forwards, vector<uint32_t> &counts) {
    for (size_t i = 0; i < path_nodes.size(); i++)
        if (forwards.at(i)) // read forward
            for(uint32_t abundance : nodes.at(path_nodes.at(i)).abundances)
                counts.push_back(abundance);
        else // read backward
            for(int k = int (nodes.at(path_nodes.at(i)).abundances.size() - 1); k > -1; k--)
                counts.push_back(nodes.at(path_nodes.at(i)).abundances.at(k));
}

bool DBG::check_path_consistency(const vector<size_t> &path_nodes, const vector<bool> &forwards) {
    // one orientation for each node
    if(path_nodes.size() != forwards.size())
        return false;

    for(size_t i = 0; i < path_nodes.size() - 1; i++){
        bool found = false;
        // is there an arc leading to a consistent node?
        for(auto &arc : nodes.at(path_nodes.at(i)).arcs)
            // same node orientation and successor check
            if(arc.forward == forwards.at(i) && arc.successor == path_nodes.at(i + 1))
                found = true;
        if(!found)
            return false;
    }
    return true;
}

uint32_t DBG::get_n_kmers() const {
    return n_kmers;
}

uint32_t DBG::get_n_nodes() const {
    return nodes.size();
}

const node_t & DBG::get_node(size_t node){
    return nodes.at(node);
}

uint32_t DBG::get_kmer_size() {
    return kmer_size;
}

