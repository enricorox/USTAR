//
// Created by enrico on 28/12/22.
//

#include <fstream>
#include <iostream>
#include "Decoder.h"
#include "DBG.h"
#include "consts.h"

Decoder::Decoder(const string &fasta_file_name, const string &counts_file_name, int kmer_size, bool debug) {
    this->fasta_file_name = fasta_file_name;
    this->counts_file_name = counts_file_name;

    this->kmer_size = kmer_size;
    this->debug = debug;
}

void Decoder::extract_kmers_and_counts(const string &output_file_name) {
    ifstream simplitigs(fasta_file_name);
    ofstream output(output_file_name);

    if(!simplitigs.good()){
        cerr << "Cannot open " << fasta_file_name << endl;
        exit(EXIT_FAILURE);
    }

    if(!output.good()){
        cerr << "Cannot open " << output_file_name << endl;
        exit(EXIT_FAILURE);
    }

    string simplitig;
    size_t n_kmers = 0;
    while(true){
        // break if end of file
        if(!getline(simplitigs, simplitig))
            break;

        // skip def-lines
        if(simplitig[0] == '>')
            continue;

        for(size_t i = 0; i < simplitig.size() - kmer_size + 1; i++){
            string kmer = simplitig.substr(i, kmer_size);
            string kmer_rc = DBG::reverse_complement(kmer);
            // output canonical kmer
            // NOTE that DSK computes canonical kmers with A<C<T<G order! Then I use jellyfish for validation
            output << ((kmer < kmer_rc) ? kmer : kmer_rc);
            output << " " << counts.at(n_kmers++) << "\n";

            if(n_kmers > counts.size()){
                cerr << "There are too few counts" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    if(n_kmers < counts.size()){
        cerr << "There are too much counts" << endl;
        exit(EXIT_FAILURE);
    }
    simplitigs.close();
    output.close();
}

void Decoder::decode(encoding_t encoding) {
    ifstream counts_file(counts_file_name);

    if(!counts_file.good()){
        cerr << "decode(): cannot open " << counts_file_name << endl;
        exit(EXIT_FAILURE);
    }

    switch(encoding){
        case encoding_t::PLAIN: {
                uint32_t c;
                while (counts_file >> c){
                    counts.push_back(c);
                }
            }
            break;
        default:
            cerr << "decode(): unknown encoding" << endl;
            exit(EXIT_FAILURE);
    }

    counts_file.close();
}
