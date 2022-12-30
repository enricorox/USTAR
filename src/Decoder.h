//
// Created by enrico on 28/12/22.
//

#ifndef USTAR_DECODER_H
#define USTAR_DECODER_H

#include <string>
#include <vector>

#include "consts.h"

using namespace std;
class Decoder{
    string fasta_file_name;
    string counts_file_name;
    int kmer_size;
    bool debug;

    vector<uint32_t> counts;

public:
    Decoder(const string &fasta_file_name, const string &counts_file_name, int kmer_size, bool debug);

    void extract_kmers_and_counts(const string &output_file_name);

    void decode(encoding_t encoding);
};
#endif //USTAR_DECODER_H
