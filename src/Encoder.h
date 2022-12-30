//
// Created by enrico on 22/12/22.
//

#ifndef USTAR_ENCODER_H
#define USTAR_ENCODER_H

#include <vector>
#include <map>
#include "consts.h"
using namespace std;

class Encoder{
    bool debug = true;
    encoding_t encoding{};
    bool encoding_done = false;

    vector<bool> flips;
    vector<size_t> simplitigs_order;
    vector<double> avg_counts;

    const vector<string> *simplitigs;
    const vector<vector<uint32_t>> *simplitigs_counts;
    size_t n_kmers = 0;

    vector<uint32_t> symbols;
    vector<uint32_t> runs;
    double avg_run = 0;

    vector<uint32_t> compacted_counts;
    vector<uint32_t> bwt_counts;
    size_t bwt_primary_index = 0;

    void do_RLE();

    void compute_avg();

    void do_flip();

    void compact_counts();

public:
    Encoder(const vector<string> *simplitigs, const vector<vector<uint32_t>> *simplitigs_counts, bool debug=false);

    void encode(encoding_t encoding_type);

    void to_fasta_file(const string &file_name);

    void to_counts_file(const string &file_name);

    void print_stat();
};
#endif //USTAR_ENCODER_H
