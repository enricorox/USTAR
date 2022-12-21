#include <iostream>
#include <getopt.h>
#include "commons.h"
#include "DBG.h"
#include "SPSS.h"

using namespace std;

void print_params(const params_t &params){
    cout << "Params:" << endl;
    cout << "\tInput file: " << params.input_file << endl;
    cout << "\tkmer_size: " << params.kmer_size << endl;
    cout << "\tverify_input: " << (params.verify_input?"true":"false") << endl;
    cout << endl;
}

void parse_cli(int argc, char **argv, params_t &params){
    int c;
    while((c = getopt(argc, argv, "i:k:v")) != -1){
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
                if(!optarg || params.kmer_size <= 0) {
                    cerr << "Need a positive kmer size!" << endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'v':
                params.verify_input = true;
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

    if(params.verify_input) {
        if (dbg.verify_overlaps())
            cout << "DBG is an overlapping graph!" << endl;
        else
            cout << "DBG is NOT an overlapping graph" << endl;
        if(dbg.validate())
            cout << "DBG is the same as BCALM2 one!" << endl;
        else
            cout << "DBG is NOT the same as BCALM2 one!" << endl;
    }

    SPSS spss(&dbg);
    spss.simpler_test();

    return EXIT_SUCCESS;
}
