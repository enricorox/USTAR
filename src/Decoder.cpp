//
// Created by enrico on 28/12/22.
//

#include <fstream>
#include <iostream>
#include "Decoder.h"
#include "DBG.h"
#include "consts.h"

Decoder::Decoder(const string &fasta_file_name, const string &counts_file_name, const string &output_file_name, int kmer_size, bool debug) {
    this->fasta_file_name = fasta_file_name;
    this->counts_file_name = counts_file_name;
    this->output_file_name = output_file_name;
    this->kmer_size = kmer_size;
    this->debug = debug;
}

void Decoder::decode() {
    ifstream simplitigs(fasta_file_name);
    ifstream counts(counts_file_name);
    ofstream output(output_file_name);

    if(!simplitigs.good()){
        cerr << "Cannot open " << fasta_file_name << endl;
        exit(EXIT_FAILURE);
    }

    if(!counts.good()){
        cerr << "Cannot open " << counts_file_name << endl;
        exit(EXIT_FAILURE);
    }

    if(!output.good()){
        cerr << "Cannot open " << output_file_name << endl;
        exit(EXIT_FAILURE);
    }

    string simplitig;
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
            int count;
            counts >> count;
            output << " " << count << "\n";
        }
    }

    simplitigs.close();
    counts.close();
    output.close();
}
