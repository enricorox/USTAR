//
// Created by enrico on 28/12/22.
//

#include <fstream>
#include <iostream>
#include "Decoder.h"
#include "DBG.h"
#include "consts.h"
#include "bwt.hpp"

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

    // line read in fasta file
    string simplitig;
    // number of kmers incremented at each written line
    size_t n_kmers = 0;
    while(true){
        // break if end of file
        if(!getline(simplitigs, simplitig))
            break;

        // skip def-lines
        if(simplitig[0] == '>')
            continue;

        // extract kmers in a simplitig
        for(size_t i = 0; i < simplitig.size() - kmer_size + 1; i++){
            string kmer = simplitig.substr(i, kmer_size);
            string kmer_rc = DBG::reverse_complement(kmer);

            // output canonical kmer only
            // NOTE that DSK computes canonical kmers with A<C<T<G order! Then I use jellyfish for validation
            output << ((kmer < kmer_rc) ? kmer : kmer_rc);

            if(n_kmers >= counts.size()){
                cerr << "extract_kmers_and_counts(): There are too few counts" << endl;
                exit(EXIT_FAILURE);
            }

            output << " " << counts.at(n_kmers++) << "\n";
        }
    }

    if(n_kmers < counts.size()){
        cerr << "extract_kmers_and_counts(): There are too much counts" << endl;
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
                // simply read integers
                uint32_t c;
                while (counts_file >> c){
                    counts.push_back(c);
                }
            }
            break;
        case encoding_t::RLE:{
                string line;
                while(getline(counts_file, line)){
                    // is this line a "S:R"?
                    size_t pos = line.find(RLE_SEPARATOR);
                    if(pos == string::npos){ // no parse_file
                        int symbol = atoi(line.c_str());
                        counts.push_back(symbol);

                        if(debug && symbol < 2){
                            cout << "Warning: found count = " << symbol << "\n";
                        }
                    } else{ // there is a parse_file here
                        // break line in two string where there is the separator
                        line[pos] = '\0';

                        // read symbol and parse_file length
                        int symbol = atoi(line.c_str());
                        int run = atoi(line.c_str() + pos + 1);

                        if(debug && run < 2)
                            cout << "Warning: found " << run << "-length parse_file \n";

                        // save multiple symbols
                        for(int i = 0; i < run; i++)
                            counts.push_back(symbol);
                    }
                }
            }
            break;
        case encoding_t::BWT: {
                counts.clear();
                // primary index is the first number
                long primary_key;
                counts_file >> primary_key;

                if(debug)
                    cout << "Primary key = " << primary_key << endl;

                // read all the other counts
                uint32_t c;
                while(counts_file >> c)
                    counts.push_back(c);

                if(debug)
                    cout << "Read " << counts.size() << " values" << endl;

                // do inverse BWT
                auto ikey = counts.begin() + primary_key;
                townsend::algorithm::bwtDecode(counts.begin(), counts.end(), ikey);
            }
            break;
        default:
            cerr << "decode(): unknown encoding" << endl;
            exit(EXIT_FAILURE);
    }
    if(debug)
        cout <<"Found " << counts.size() << " counts\n";
    counts_file.close();
}
