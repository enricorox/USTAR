//
// Created by enrico on 24/01/23.
//

#include <fstream>
#include <iostream>
#include "Analyzer.h"

Analyzer::Analyzer(const string &file_name, int kmer_size) {
    this->file_name = file_name;
    this->kmer_size = kmer_size;

    parse_file();
    compute_variance();
}

void Analyzer::parse_file() {
    ifstream in(file_name);

    if(!in.good()){
        cerr << "parse_file(): cannot open " << file_name << endl;
        exit(EXIT_FAILURE);
    }

    F[A] = F[C] = F[T] = F[G] = 0;
    string line;
    while(getline(in, line)){
        int i;
        for(i = 0; i < kmer_size; i++){
            switch(line[i]){
                case 'A':
                case 'a':
                    F[A]++;
                    break;
                case 'C':
                case 'c':
                    F[C]++;
                    break;
                case 'T':
                case 't':
                    F[T]++;
                    break;
                case 'G':
                case 'g':
                    F[G]++;
                    break;
                default:
                    cerr << "parse_file(): unexpected character: " << line[i] << endl;
                    exit(EXIT_FAILURE);
            }
        }

        if(line[i] != ' '){
            cerr << "parse_file(): blank expected at column " << i << endl;
            exit(EXIT_FAILURE);
        }

        int ab;
        sscanf(line.c_str() + i + 1, "%d", &ab);
        counts.push_back(ab);
    }

}

void Analyzer::print_stats() {
    cout << "kmers stats:\n";
    cout << "   variance:   " << variance << "\n";
    cout << "   GC-content:" << double (F[C] + F[G]) / (F[A] + F[C] + F[T] + F[G]) << "\n";
    cout << endl;
}

void Analyzer::compute_variance() {
    // first compute mean
    double sum = 0;
    for(auto c : counts)
        sum += c;
    double mean = sum / double (counts.size());

    // then compute variance
    double variance_sum = 0;
    for(auto c : counts)
        variance_sum += (c - mean) * (c - mean);

    variance = variance_sum / double (counts.size());
}
