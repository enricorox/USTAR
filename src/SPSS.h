//
// Created by enrico on 21/12/22.
//

#ifndef USTAR_SPSS_H
#define USTAR_SPSS_H

#include "DBG.h"

using namespace std;

class SPSS{
    DBG *dbg;
    uint32_t n_nodes;

    // simplitigs
    vector<vector<size_t>> simplitigs_path_nodes;
    vector<vector<bool>> simplitigs_path_forwards;
    size_t n_simplitigs = 0;

    // saturated nodes
    vector<bool> saturated;

    void extends(uint32_t seed, vector<size_t> &path_nodes, vector<bool> &path_forwards);

public:
    explicit SPSS(DBG *dbg);

    void extract_simplitigs();

    void to_fasta_file(const string &file_name);

    void to_counts_file(const string &file_name);

    void test();
};

#endif //USTAR_SPSS_H
