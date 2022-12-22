#include <iostream>
#include <getopt.h>
#include "commons.h"
#include "DBG.h"
#include "SPSS.h"

using namespace std;

void print_help(){
    cout << "Usage: ./USTAR -i <bcalm_file> -k <kmer_size>" << endl;
    cout << "Options:" << endl;
    cout << "\t-o" << "\toutput file" << endl;
    cout << "\t-v" << "\tvalidate input" << endl;
}

void print_params(const params_t &params){
    cout << "Params:" << endl;
    cout << "\tinput file: " << params.input_file << endl;
    cout << "\toutput file: " << params.output_file << endl;
    cout << "\tkmer_size: " << params.kmer_size << endl;
    cout << "\tverify_input: " << (params.verify_input?"true":"false") << endl;
    cout << endl;
}

void parse_cli(int argc, char **argv, params_t &params){
    bool done = false;
    int c;
    while((c = getopt(argc, argv, "i:k:vo:")) != -1){
        done = true;
        switch(c){
            case 'i':
                if(optarg)
                    params.input_file = string(optarg);
                else {
                    cerr << "Need an input file name!" << endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'o':
                if(optarg)
                    params.output_file = string(optarg);
                else {
                    cerr << "Need an output file name!" << endl;
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
                print_help();
                exit(EXIT_FAILURE);
        }
    }
    if(!done){
        print_help();
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char **argv) {
    cout << "===== USTAR by enricorox =====" << endl;
    // cli parameters
    params_t params;
    parse_cli(argc, argv, params);
    print_params(params);

    // make a dBG
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

    // make an SPSS
    SPSS spss(&dbg);

    // compute simplitigs
    spss.extract_simplitigs();
    spss.to_fasta_file(params.output_file + ".ustar.fa");
    spss.to_counts_file(params.output_file + ".ustar.counts");

    return EXIT_SUCCESS;
}
