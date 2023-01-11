#include <iostream>
#include <getopt.h>
#include <chrono>

#include "DBG.h"
#include "SPSS.h"
#include "Encoder.h"
#include "consts.h"
#include "commons.h"

using namespace std;
using namespace std::chrono;

struct params_t{
    string input_file_name{};
    string output_file_name = "out";
    string fasta_file_name = output_file_name + ".ustar.fa";
    string counts_file_name = output_file_name + ".ustar.counts";

    int kmer_size = 31;

    bool debug = false;

    encoding_t encoding = encoding_t::PLAIN;
    seeding_method_t seeding_method = seeding_method_t::FIRST;
    extending_method_t extending_method = extending_method_t::FIRST;
    int iterations = 1;
};


void print_help(const params_t &params){
    cout << "Find a Spectrum Preserving String Set (aka simplitigs) for the input file.\n";
    cout << "Compute the kmer counts vector.\n\n";

    cout << "Usage: ./USTAR -i <input_file_name>\n\n";

    cout << "Basic options:\n\n";

    cout << "   -k  kmer size, must be the same of BCALM2 [" << params.kmer_size << "]\n\n";

    cout << "   -c  counts file name [" << params.counts_file_name << "]\n\n";

    cout << "   -o  output file base name [" << params.output_file_name << "]\n\n";

    cout << "   -v  print version and author\n\n";

    cout << "   -h  print this help\n\n" << endl;

    cout << "Advanced options:\n\n";

    cout << "   -s  seeding method [" << inv_map<seeding_method_t>(seeding_method_names, params.seeding_method) << "]\n";
    cout << "       f               choose the first seed available\n";
    cout << "       r               choose a random seed\n";
    cout << "       -ma             choose the seed with lower median abundance\n";
    cout << "       -aa             choose the seed with lower average abundance\n";
    cout << "       =a              choose the seed with most similar abundance to the last seòected node\n";
    cout << "       -l              choose the seed with smaller length\n";
    cout << "       +l              choose the seed with bigger length\n";
    cout << "       -c              choose the seed with less arcs\n";
    cout << "       +c              choose the seed with more arcs\n";
    cout << "\n";

    cout << "   -x  extending method [" << inv_map<extending_method_t>(extending_method_names, params.extending_method) << "]\n";
    cout << "       f               choose the first successor available\n";
    cout << "       r               choose a random successor\n";
    cout << "       =a              choose the successor with most similar abundance to the last selected node\n";
    cout << "       -l              choose the successor with smaller length\n";
    cout << "       +l              choose the successor with bigger length\n";
    cout << "       -c              choose the successor with less arcs\n";
    cout << "       +c              choose the successor with more arcs\n";
    cout << "\n";

    cout << "   -I  number of iterations [" << params.iterations << "]\n";
    cout << "       available only for random methods\n\n";

    cout << "   -e  encoding [" << inv_map<encoding_t>(encoding_names, params.encoding)<< "]\n";
    cout << "       plain           do not use any encoding\n";
    cout << "       rle             use special Run Length Encoding\n";
    cout << "       avg_rle         sort simplitigs by average counts and use RLE\n";
    cout << "       flip_rle        make contiguous runs by flipping simplitigs if necessary and use RLE\n";
    cout << "       avg_flip_rle    make contiguous runs by sorting by average, flipping simplitigs if necessary and use RLE\n";
    cout << "\n";

    cout << "   -d  debug [" << (params.debug?"true":"false") << "]\n\n";
}

void print_params(const params_t &params){
    cout << "Params:\n";
    cout << "   input file:             " << params.input_file_name << "\n";
    cout << "   kmer size:              " << params.kmer_size << "\n";
    cout << "   output file base name:  " << params.output_file_name << "\n";
    cout << "   counts file name:       " << params.counts_file_name << "\n";
    cout << "   seeding method:         " << inv_map<seeding_method_t>(seeding_method_names, params.seeding_method) << "\n";
    cout << "   extending method:       " << inv_map<extending_method_t>(extending_method_names, params.extending_method) << "\n";
    cout << "   iterations:             " << params.iterations << "\n";
    cout << "   encoding:               " << inv_map<encoding_t>(encoding_names, params.encoding) << "\n";
    cout << "   debug:                  " << (params.debug?"true":"false") << "\n";
    cout << endl;
}

void parse_cli(int argc, char **argv, params_t &params){
    bool done = false;
    int c;
    while((c = getopt(argc, argv, "i:k:vo:dhe:s:x:c:I:")) != -1){
        switch(c){
            case 'i':
                params.input_file_name = string(optarg);
                done = true;
                break;
            case 'o':
                params.output_file_name = string(optarg);
                params.fasta_file_name = params.output_file_name + ".ustar.fa";
                params.counts_file_name =
                        params.output_file_name + ".ustar" + encoding_suffixes.at(params.encoding) + ".counts";
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
                params.counts_file_name = params.output_file_name
                                          + ".ustar" + encoding_suffixes.at(params.encoding) + ".counts";
                break;
            case 's': // seed method
                if(seeding_method_names.find(optarg) == seeding_method_names.end()){
                    cerr << "parse_cli(): " << optarg << " is not a valid seed method" <<endl;
                    exit(EXIT_FAILURE);
                }
                params.seeding_method = seeding_method_names.at(optarg);
                break;
            case 'x': // extension method
                if(extending_method_names.find(optarg) == extending_method_names.end()){
                    cerr << "parse_cli(): " << optarg << " is not a valid extension method" <<endl;
                    exit(EXIT_FAILURE);
                }
                params.extending_method = extending_method_names.at(optarg);
                break;
            case 'I':
                params.iterations = atoi(optarg);
                if(params.iterations <= 0) {
                    cerr << "parse_cli(): Need a positive number of iterations!" << endl;
                    exit(EXIT_FAILURE);
                }
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

int main(int argc, char **argv) {
    cout << "===== Unitig STitch STar (USTAR) v" << VERSION << " =====\n";
    // cli parameters
    params_t params;
    parse_cli(argc, argv, params);
    print_params(params);

    // make a dBG
    cout << "Reading the input file..." << endl;
    auto start_time = steady_clock::now();
    DBG dbg(params.input_file_name, params.kmer_size, params.debug);
    auto stop_time = steady_clock::now();
    cout << "Reading time: " << duration_cast<seconds>(stop_time - start_time).count() << " s\n";
    dbg.print_stat();

    // verify input
    if(params.debug)
        dbg.verify_input();


    // choose SPSS sorter
    Sorter sorter(params.seeding_method, params.extending_method, params.debug);
    // make an SPSS
    SPSS spss(&dbg, &sorter, params.debug);

    cout << "Computing a path cover..." << endl;
    start_time = steady_clock::now();
    spss.compute_path_cover();
    stop_time = steady_clock::now();
    cout << "Computing time: " << duration_cast<milliseconds>(stop_time - start_time).count() << " ms\n";

    cout << "Extracting simplitigs and kmers counts..." << endl;
    spss.extract_simplitigs_and_counts();
    spss.print_stats();

    Encoder encoder(spss.get_simplitigs(), spss.get_counts(), params.debug);
    encoder.encode(params.encoding);
    encoder.print_stat();
    encoder.to_fasta_file(params.fasta_file_name);
    cout << "Simplitigs written to disk: " << params.fasta_file_name << endl;
    encoder.to_counts_file(params.counts_file_name);
    cout << "Counts written to disk: " << params.counts_file_name << endl;

    return EXIT_SUCCESS;
}
