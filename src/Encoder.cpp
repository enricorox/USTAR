//
// Created by enrico on 22/12/22.
//

#include <fstream>
#include <iostream>
#include <algorithm>
#include "Encoder.h"

Encoder::Encoder(const vector<string> *simplitigs, const vector<vector<uint32_t>> *counts, bool debug) {
    this->debug = debug;
    this->simplitigs = simplitigs;
    this->counts = counts;

    simplitigs_order.resize(simplitigs->size());
    for(size_t i = 0; i < simplitigs->size(); i++)
        simplitigs_order.push_back(i);
}

void Encoder::to_fasta_file(const string &file_name) {
    if(simplitigs->empty()){
        cerr << "to_fasta_file(): There are no simplitigs!" << endl;
        exit(EXIT_FAILURE);
    }

    ofstream fasta;
    fasta.open(file_name);
    for(size_t i = 0; i < simplitigs->size(); i++){
        auto &simplitig = (*simplitigs)[simplitigs_order[i]];
        fasta << ">\n";
        fasta << simplitig << "\n";
    }
    fasta.close();
}

void Encoder::to_counts_file(const string &file_name) {
    if (counts->empty()) {
        cerr << "to_counts_file(): There are no counts!" << endl;
        exit(EXIT_FAILURE);
    }

    ofstream encoded;
    switch(encoding) {
        case AVG_RLE:
            // no break here
        case RLE:
            encoded.open(file_name + ".rle");
            if(!encoded.good()){
                cerr << "Can't open output file: " << file_name + ".rle" << endl;
                exit(EXIT_FAILURE);
            }

            for(size_t i = 0; i < symbols.size(); i++){
                encoded << symbols[i];
                if(runs[i] != 1)
                    encoded << " " << runs[i];
                encoded << "\n";
            }
            encoded.close();
            break;

        case PLAIN:
        default:
            encoded.open(file_name);
            if(!encoded.good()){
                cerr << "Can't open output file:" << file_name << endl;
                exit(EXIT_FAILURE);
            }
            for (const auto &simplitig_counts: *counts) {
                for (auto c: simplitig_counts)
                    encoded << c << (debug ? " " : "\n");
                if (debug) encoded << "\n";
            }
    }
    encoded.close();
}

void Encoder::encode(encoding_t encoding_type) {
    this->encoding = encoding_type;

    switch(encoding) {
        case RLE:
            do_RLE();
            break;
        case AVG_RLE:
            compute_avg();
            sort(simplitigs_order.begin(), simplitigs_order.end(),
                 [this](size_t a, size_t b){return avg_counts[a] < avg_counts[b];}
                 );
            do_RLE();
            break;
        case PLAIN:
        default:
            break;
    }
}

void Encoder::do_RLE(){
    uint32_t c = 0, prev;
    bool first = true;

    for (size_t i = 0; i < counts->size(); i++) {
        for (uint32_t curr: (*counts)[simplitigs_order[i]]) {
            // stream of counts here
            if (first || curr == prev) {
                c++;
                first = false;
            } else {
                symbols.push_back(prev);
                runs.push_back(c);
                // reset
                c = 1;
            }
            prev = curr;
        }
    }
    // save last run
    symbols.push_back(prev);
    runs.push_back(c);
}

void Encoder::compute_avg() {
    avg_counts.resize(counts->size());

    for(auto &simplitig_counts : *counts){
        double sum = 0;
        for(uint32_t c : simplitig_counts)
            sum += c;
        avg_counts.push_back(sum / (double) simplitig_counts.size());
    }
}

