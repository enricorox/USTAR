//
// Created by enrico on 22/12/22.
//

#include <fstream>
#include <iostream>
#include "Encoder.h"

Encoder::Encoder(const vector<string> *simplitigs, const vector<vector<uint32_t>> *counts) {
    this->simplitigs = simplitigs;
    this->counts = counts;
}

void Encoder::to_fasta_file(const string &file_name) {
    if(simplitigs->empty()){
        cerr << "to_fasta_file(): There are no simplitigs!" << endl;
        exit(EXIT_FAILURE);
    }

    ofstream fasta;
    fasta.open(file_name);
    for(auto &simplitig : *simplitigs){
        fasta << ">\n";
        fasta << simplitig << "\n";
    }
    fasta.close();
}

void Encoder::to_counts_file(const string &file_name) {
    if(counts->empty()){
        cerr << "to_counts_file(): There are no counts!" << endl;
        exit(EXIT_FAILURE);
    }

    ofstream counts_file;
    counts_file.open(file_name);
    for(const auto &simplitig_counts : *counts){
        for(auto c : simplitig_counts)
            counts_file << c << (debug?" ":"\n");
        if(debug) counts_file << "\n";
    }
    counts_file.close();
}


