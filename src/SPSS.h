//
// Created by enrico on 21/12/22.
//

#ifndef USTAR_SPSS_H
#define USTAR_SPSS_H

#include "DBG.h"
#include "Encoder.h"

using namespace std;

class SPSS{
    bool debug;

    DBG *dbg;
    uint32_t n_nodes;

    // simplitigs
    vector<vector<size_t>> simplitigs_path_nodes;
    vector<vector<bool>> simplitigs_path_forwards;
    size_t n_simplitigs = 0;

    // saturated nodes
    vector<bool> saturated;

    void extends(uint32_t seed, vector<size_t> &path_nodes, vector<bool> &path_forwards, bool two_way);

public:
    explicit SPSS(DBG *dbg, bool debug=false);

    void extract_simplitigs(bool two_way=true);

    void to_fasta_file(const string &file_name);

    void to_counts_file(const string &file_name);

    void print_info();
};

#endif //USTAR_SPSS_H
