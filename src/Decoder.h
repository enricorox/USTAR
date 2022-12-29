//
// Created by enrico on 28/12/22.
//

#ifndef USTAR_DECODER_H
#define USTAR_DECODER_H

#include <string>

using namespace std;
class Decoder{
    string fasta_file_name;
    string counts_file_name;
    string output_file_name;
    int kmer_size;
    bool debug;

public:
    Decoder(string fasta_file_name, string counts_file_name, string output_file_name, int kmer_size, bool debug);

    void decode();
};
#endif //USTAR_DECODER_H
