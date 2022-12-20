#include <iostream>
#include <getopt.h>
#include "commons.h"
#include "DBG.h"

using namespace std;

void print_params(const params_t &params){
    cout << "Input file: " << params.input_file << endl;
    cout << "kmer-size: " << params.kmer_size << endl;
}

void parse_cli(int argc, char **argv, params_t &params){
    int c;
    while((c = getopt(argc, argv, "i:k:")) != -1){
        switch(c){
            case 'i':
                if(optarg)
                    params.input_file = string(optarg);
                else {
                    cerr << "Need a file name!" << endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'k':
                if(optarg)
                    params.kmer_size = atoi(optarg);
                else {
                    cerr << "Need a file name!" << endl;
                    exit(EXIT_FAILURE);
                }
                break;
            default:
                cerr << "WARNING: unknown parameter!" << endl;
        }
    }
}

int main(int argc, char **argv) {
    cout << "===== USTAR by enricorox =====" << endl;
    params_t params;
    parse_cli(argc, argv, params);
    print_params(params);

    DBG dbg(params.input_file, params.kmer_size);
    dbg.print_info();
    dbg.verify_overlaps();
    return 0;
}
