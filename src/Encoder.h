//
// Created by enrico on 22/12/22.
//

#ifndef USTAR_ENCODER_H
#define USTAR_ENCODER_H

#include <vector>

using namespace std;

class Encoder{
    const vector<string> *simplitigs;
    const vector<vector<size_t>> *counts;
    bool debug = true;

public:
    Encoder(const vector<string> *simplitigs, const vector<vector<size_t>> *counts);

    void to_fasta_file(const string &file_name);

    void to_counts_file(const string &file_name);
};
#endif //USTAR_ENCODER_H
