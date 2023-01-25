//
// Created by enrico on 24/01/23.
//

#ifndef USTAR_ANALYZER_H
#define USTAR_ANALYZER_H

#include <string>
#include <vector>

using namespace std;

class Analyzer{
    const int A = 0, C = 1, T = 2, G = 3;
    long F[4];

    double variance;

    vector<int> counts;

    string file_name;
    int kmer_size;

    void parse_file();

    void compute_variance();

public:
    Analyzer(const string &file_name, int kmer_size);

    void print_stats();
};


#endif //USTAR_ANALYZER_H
