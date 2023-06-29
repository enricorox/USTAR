//
// Created by enrico on 21/12/22.
//

#ifndef USTAR_SPSS_H
#define USTAR_SPSS_H

#include "DBG.h"
#include "Sorter.h"
#include "Encoder.h"

using namespace std;

class SPSS{
    bool debug;
    bool duplicates;

    DBG *dbg;
    size_t n_nodes;

    Sorter *sorter;

    // path cover
    vector<vector<node_idx_t>> path_cover_nodes;
    vector<vector<bool>> path_cover_forwards;
    // visited nodes
    vector<bool> visited;

    // simplitigs
    vector<string> simplitigs;
    vector<vector<uint32_t>> counts;
    size_t n_simplitigs = 0;

    void extends(vector<node_idx_t> &path_nodes, vector<bool> &path_forwards);

public:
    SPSS(DBG *dbg, Sorter *sorter, bool duplicates=false, bool debug=false);

    void compute_path_cover(bool two_way=true);

    void print_stats();

    void extract_simplitigs_and_counts();

    const vector<string> * get_simplitigs();

    const vector<vector<uint32_t>> * get_counts();

    size_t get_score();
};

#endif //USTAR_SPSS_H
