//
// Created by enrico on 24/01/23.
//

#include <iostream>
#include <getopt.h>
#include "consts.h"
#include "Analyzer.h"

using namespace std;

struct params_t{
    string file_name{};
    int kmer_size = 31;
};

void print_params(const params_t &params){
    cout << "Params:\n";
    cout << "   file name:  " << params.file_name << "\n";
    cout << "   kmer size:  " << params.kmer_size << "\n";
    cout << endl;
}

void print_help(const params_t &params){
    cout << "Print kmers file stats\n";
    cout << "Usage: ./ustars -i <kmers-file> -k <kmer-size>\n";
    cout << endl;
}

void parse_cli(int argc, char **argv, params_t &params){
    bool got_input = false;
    int c;
    while((c = getopt(argc, argv, "i:k:hv")) != -1) {
        switch (c) {
            case 'i':
                params.file_name = string(optarg);
                got_input = true;
                break;
            case 'k':
                params.kmer_size = atoi(optarg);
                break;
            case 'h':
                print_help(params);
                exit(EXIT_SUCCESS);
            case '?':
                cerr << "parse_cli(): missing argument or invalid option\n\n";
                print_help(params);
                exit(EXIT_FAILURE);
            default: // should never go here
                cerr << "parse_cli(): unknown option in optstring '" << c << "'\n\n";
                print_help(params);
                exit(EXIT_FAILURE);
        }
    }

    if(!got_input){
        cerr << "parse_cli(): mandatory parameter -i" << endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char **argv){
    cout << "===== Unitig STitch STar Statistics (USTARS) v" << VERSION << " =====\n";

    params_t params;
    parse_cli(argc, argv, params);
    print_params(params);
    cout << "Reading file...\n";
    Analyzer analyzer(params.file_name, params.kmer_size);
    analyzer.print_stats();
}