//
// Created by enrico on 20/12/22.
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include "DBG.h"

void DBG::parse_bcalm_file(const string &bcalm_file_name) {
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
    parse_bcalm_file(bcalm_file_name);

    for(auto &node : nodes) {
        n_edges += node.edges.size();
        n_kmers += node.abundances.size();
    }
}

DBG::~DBG() {

}

void DBG::print_info() {
    cout << "Info for " << bcalm_file_name << ":" << endl;
    cout << "\tnumber of kmers: " << n_kmers << endl;
    cout << "\tnumber of nodes: " << nodes.size() << endl;
    cout << "\tnumber of edges: " << n_edges << endl;
    cout << "\taverage number of edges: " << 1.0 * n_edges / nodes.size() << endl;
}

bool DBG::verify_overlaps() {
    for(auto &node : nodes){
        for(auto &edge : node.edges)
            if(!overlaps(node, edge))
                return false;
    }
    return true;
}

bool DBG::overlaps(const node_t &node, const edge_t &edge){

}