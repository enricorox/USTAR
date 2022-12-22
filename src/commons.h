//
// Created by enrico on 20/12/22.
//

#ifndef USTAR_COMMONS_H
#define USTAR_COMMONS_H

using namespace std;

struct params_t{
    string input_file = "../experiments/k31.a1.unitigs.fa";
    string output_file = "out";
    int kmer_size = 31;
    bool debug = false;
};

#endif //USTAR_COMMONS_H
