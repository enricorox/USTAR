//
// Created by enrico on 22/12/22.
//

#include <fstream>
#include <iostream>
#include <algorithm>
#include "Encoder.h"

Encoder::Encoder(const vector<string> *simplitigs, const vector<vector<uint32_t>> *simplitigs_counts, bool debug) {
    if (simplitigs_counts->empty()) {
        cerr << "to_counts_file(): There are no simplitigs_counts!" << endl;
        exit(EXIT_FAILURE);
    }
    if(simplitigs->size() != simplitigs_counts->size()){
        cerr << "Encoder(): Got an incorrect number of simplitigs_counts or simplitigs!" << endl;
        exit(EXIT_FAILURE);
    }

    this->debug = debug;
    this->simplitigs = simplitigs;
    this->simplitigs_counts = simplitigs_counts;

    simplitigs_order.reserve(simplitigs->size());
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
    ofstream encoded;
    encoded.open(file_name);
    if(!encoded.good()){
        cerr << "Can't open output file: " << file_name + ".rle" << endl;
        exit(EXIT_FAILURE);
    }

    switch(encoding) {
        case AVG_RLE:
            // no break here
        case RLE:
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
            for (const auto &simplitig_counts: *simplitigs_counts) {
                for (auto c: simplitig_counts)
                    encoded << c << (debug ? " " : "\n");
                if (debug) encoded << "\n";
            }
    }
    encoded.close();
}

void Encoder::encode(encoding_t encoding_type) {
    this->encoding = encoding_type;
    encoding_done = true;

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
    uint32_t count = 0, prev;
    uint32_t sum_run = 0;

    bool first = true;
    for (size_t i = 0; i < simplitigs_counts->size(); i++) {
        // accessing simplitigs_counts based on computed order!
        for (uint32_t curr: (*simplitigs_counts)[simplitigs_order[i]]) {
            // stream of simplitigs_counts here
            if (first || curr == prev) {
                count++;
                first = false;
            } else {
                symbols.push_back(prev);
                runs.push_back(count);
                sum_run += count;
                // reset
                count = 1;
            }
            prev = curr;
        }
    }
    // save last run
    symbols.push_back(prev);
    runs.push_back(count);
    sum_run += count;

    cout << "sum_run = " << sum_run << endl;
    cout << "run.size() = " << runs.size() << endl;

    avg_run = (double) sum_run / (double) runs.size();
}

void Encoder::compute_avg() {
    // pre-allocate the vector
    avg_counts.reserve(simplitigs_counts->size());

    for(auto &simplitig_counts : *simplitigs_counts){
        double sum = 0;
        for(uint32_t c : simplitig_counts)
            sum += c;
        avg_counts.push_back(sum / (double) simplitig_counts.size());
    }
}

void Encoder::print_stat(){
    if(!encoding_done){
        cerr << "print_stat(): Need to encode() first!" << endl;
        exit(EXIT_FAILURE);
    }
    cout << "\nEncoding stats:\n";
    switch (encoding) {
        case AVG_RLE:
            // no break here
        case RLE:
            cout << "\tNumber of runs: " << runs.size() << "\n";
            cout << "\tAverage run: " << avg_run << "\n";
            break;
        case PLAIN:
            // no break here
        default:
            cout << "\tNo size reduction for simplitigs_counts\n";
            break;
    }
    cout << endl;
}
