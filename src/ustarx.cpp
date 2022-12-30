//
// Created by enrico on 29/12/22.
//

#include <cstdlib>
#include <iostream>
#include <getopt.h>
#include <chrono>

#include "commons.h"
#include "consts.h"
#include "Decoder.h"

using namespace std;
using namespace std::chrono;

struct params_t{
    string input_file_name = "out.ustar.fa";
    string counts_file_name = "out.ustar.counts";
    string output_file_name = "ustar-kmers.txt";
    encoding_t encoding = encoding_t::PLAIN;
    int kmer_size = 31;
    bool sort = false;
    bool debug = false;
};

void print_help(const params_t &params){
    cout << "Decode and extract canonical kmers and counts from USTAR (or UST) output files.\n\n";

    cout << "Syntax: ./ustarx -i <ustar-fasta> -c <ustar-counts>\n\n";

    cout << "Options:\n\n";

    cout << "   -k  kmer size [" << params.kmer_size << "]\n\n";

    cout << "   -o  output file name [" << params.output_file_name << "]\n\n";

    cout << "   -s  sort output file [" << (params.sort ? "true" : "false") << "]\n\n";

    cout << "   -e  encoding [" << inv_map<encoding_t>(encoding_names, params.encoding)<< "]\n";
    cout << "       plain           do not use any encoding\n";
    cout << "       rle             use special Run Length Encoding\n";
    cout << "\n";

    cout << "   -d  debug [" << (params.debug ? "true" : "false") << "]\n\n";

    cout << "   -v  print author and version\n\n";

    cout << "   -h  print this help\n\n";

    cout << endl;
}

void print_params(const params_t &params){
    cout << "Params:\n";
    cout << "   input file:             " << params.input_file_name << "\n";
    cout << "   counts file name:       " << params.counts_file_name << "\n";
    cout << "   kmer size:              " << params.kmer_size << "\n";
    cout << "   output file name:       " << params.output_file_name << "\n";
    cout << "   sort:                   " << (params.sort ? "true" : "false") << "\n";
    cout << "   encoding:               " << inv_map<encoding_t>(encoding_names, params.encoding) << "\n";
    cout << "   debug:                  " << (params.debug?"true":"false") << "\n";
    cout << endl;
}

void parse_cli(int argc, char **argv, params_t &params){
    bool done = false;
    int c;
    while((c = getopt(argc, argv, "i:k:vo:dhe:c:s")) != -1){
        switch(c){
            case 'i':
                params.input_file_name = string(optarg);
                done = true;
                break;
            case 'o':
                params.output_file_name = string(optarg);
                break;
            case 'c':
                params.counts_file_name = string(optarg);
                break;
            case 'k':
                params.kmer_size = atoi(optarg);
                if(params.kmer_size <= 0) {
                    cerr << "parse_cli(): Need a positive kmer size!" << endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'v':
                cout << "Version: " << VERSION << "\n";
                cout << "Author: Enrico Rossignolo <enricorrx at gmail dot com>" << endl;
                exit(EXIT_SUCCESS);
            case 'd':
                params.debug = true;
                break;
            case 'e': // encoding
                // is a valid encoding?
                if(encoding_names.find(optarg) == encoding_names.end()){
                    cerr << "parse_cli(): " << optarg << " is not a valid encoding" <<endl;
                    exit(EXIT_FAILURE);
                }
                params.encoding = encoding_names.at(optarg);
                break;
            case 's':
                params.sort = true;
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
    if(!done){
        print_help(params);
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char **argv){
    cout << "===== USTAR kmer eXtractor (USTARX) v" << VERSION << " =====\n";

    params_t params;
    parse_cli(argc, argv, params);
    print_params(params);

    Decoder decoder(params.input_file_name, params.counts_file_name, params.kmer_size, params.debug);

    cout << "Decoding counts...\n";
    auto start_time = steady_clock::now();
    decoder.decode(params.encoding);
    auto stop_time = std::chrono::steady_clock::now();
    cout << "Done. Decoding time: " << duration_cast<milliseconds>(stop_time - start_time).count() << " ms\n";

    cout << "Extracting kmers...\n";
    start_time = steady_clock::now();
    decoder.extract_kmers_and_counts(params.output_file_name);
    stop_time = std::chrono::steady_clock::now();
    cout << "Done. Extraction time: " << duration_cast<seconds>(stop_time - start_time).count() << " s\n";

    if(params.sort){
        cout << "Sorting kmers...\n";
        start_time = steady_clock::now();
        if(system(("sort " + params.output_file_name + " -o " + params.output_file_name).c_str()) == EXIT_SUCCESS) {
            stop_time = steady_clock::now();
            cout << "Done. Sorting time: " << duration_cast<seconds>(stop_time - start_time).count() << " s\n";
        }else{
            cerr << "Something went wrong while sorting\n";
            exit(EXIT_FAILURE);
        }
    }

    exit(EXIT_SUCCESS);
}