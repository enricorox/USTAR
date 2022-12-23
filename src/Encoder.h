//
// Created by enrico on 22/12/22.
//

#ifndef USTAR_ENCODER_H
#define USTAR_ENCODER_H

#include <vector>

using namespace std;

enum encoding_t{
    PLAIN,
    RLE
};

class Encoder{
    bool debug = true;
    encoding_t encoding{};

    const vector<string> *simplitigs;
    const vector<vector<uint32_t>> *counts;

    vector<uint32_t> symbols;
    vector<uint32_t> runs;

public:
    Encoder(const vector<string> *simplitigs, const vector<vector<uint32_t>> *counts, bool debug=false);

    void encode(encoding_t encoding_type);

    void to_fasta_file(const string &file_name);

    void to_counts_file(const string &file_name);
};
#endif //USTAR_ENCODER_H
