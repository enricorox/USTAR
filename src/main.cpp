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
    cout << "\t-v" << "\tvalidate input and exit" << endl;
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
    cout << "===== Unitig STitch STar (USTAR) =====" << endl;
    // cli parameters
    params_t params;
    parse_cli(argc, argv, params);
    print_params(params);

    // make a dBG
    DBG dbg(params.input_file, params.kmer_size);
    dbg.print_info();

    if(params.verify_input) {
        if(dbg.verify_input())
            exit(EXIT_SUCCESS);
        else
            exit(EXIT_FAILURE);
    }

    // make an SPSS
    SPSS spss(&dbg);

    // compute simplitigs
    cout << "Extracting simplitigs..." << endl;
    spss.extract_simplitigs();

    spss.print_info();

    // save to disk
    string fasta_file_name = params.output_file + ".ustar.fa";
    string counts_file_name = params.output_file + ".ustar.counts";
    spss.to_fasta_file(fasta_file_name);
    cout << "Simplitigs written to disk: " << fasta_file_name << endl;
    spss.to_counts_file(counts_file_name);
    cout << "Counts written to disk: " << counts_file_name << endl;

    return EXIT_SUCCESS;
}
