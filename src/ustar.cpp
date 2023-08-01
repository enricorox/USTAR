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
    string fasta_file_name{};
    string counts_file_name{};

    int kmer_size = 31;

    int visited_depth = 0;

    bool debug = false;

    encoding_t encoding = encoding_t::PLAIN;
    seeding_method_t seeding_method = seeding_method_t::FIRST;
    extending_method_t extending_method = extending_method_t::FIRST;
    bool phases = false;
    bool profile = false;
};


void print_help(const params_t &params){
    cout << "Find a Spectrum Preserving String Set (aka simplitigs) for the input file.\n";
    cout << "Compute the kmer counts vector.\n\n";

    cout << "Usage: ./USTAR -i <input_file_name>\n\n";

    cout << "Basic options:\n\n";

    cout << "   -k  kmer size, must be the same of BCALM2 [" << params.kmer_size << "]\n\n";

    cout << "   -c  counts file name [" << params.counts_file_name << "]\n\n";

    cout << "   -o  fasta file name [" << params.fasta_file_name << "]\n\n";

    cout << "   -D  depth for exploring visited nodes [" << params.visited_depth << "]\n\n";

    cout << "   -p  connect visited nodes in a second phase [" << (params.phases ? "true":"false") << "]\n\n";

    cout << "   -P  only profile dBG [" << (params.profile ? "true":"false") << "]\n\n";

    cout << "   -v  print version and author\n\n";

    cout << "   -h  print this help\n\n" << endl;

    cout << "Advanced options:\n\n";

    cout << "   -s  seeding method [" << inv_map<seeding_method_t>(seeding_method_names, params.seeding_method) << "]\n";
    cout << "       f               choose the first seed available\n";
    cout << "       r               choose a random seed\n";
    cout << "       -ma             choose the seed with lower median abundance\n";
    cout << "       +aa             choose the seed with higher average abundance\n";
    cout << "       -aa             choose the seed with lower average abundance\n";
    cout << "       =a              choose the seed with most similar abundance to the last selected node\n";
    cout << "       -l              choose the seed with smaller length\n";
    cout << "       +l              choose the seed with bigger length\n";
    cout << "       -c              choose the seed with less arcs\n";
    cout << "       +c              choose the seed with more arcs\n";
    cout << "       -u              choose the seed less unbalanced\n";
    cout << "       +u              choose the seed more unbalanced\n";
    cout << "\n";

    cout << "   -x  extending method [" << inv_map<extending_method_t>(extending_method_names, params.extending_method) << "]\n";
    cout << "       f               choose the first successor available\n";
    cout << "       r               choose a random successor\n";
    cout << "       =a              choose the successor with most similar abundance to the last selected node\n";
    cout << "       =ma             choose the successor with most similar median abundance to the last selected node\n";
    cout << "       -ma             choose the successor with lower abundance to the last selected node\n";
    cout << "       -l              choose the successor with smaller length\n";
    cout << "       +l              choose the successor with bigger length\n";
    cout << "       -c              choose the successor with less arcs\n";
    cout << "       +c              choose the successor with more arcs\n";
    cout << "\n";

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
    cout << "   visited nodes depth:    " << params.visited_depth << "\n";
    cout << "   phases:                 " << (params.phases?"true":"false") << "\n";
    cout << "   fasta file name:        " << params.fasta_file_name << "\n";
    cout << "   counts file name:       " << params.counts_file_name << "\n";
    cout << "   seeding method:         " << inv_map<seeding_method_t>(seeding_method_names, params.seeding_method) << "\n";
    cout << "   extending method:       " << inv_map<extending_method_t>(extending_method_names, params.extending_method) << "\n";
    cout << "   encoding:               " << inv_map<encoding_t>(encoding_names, params.encoding) << "\n";
    cout << "   profile:                " << (params.profile?"true":"false") << "\n";
    cout << "   debug:                  " << (params.debug?"true":"false") << "\n";
    cout << endl;
}

void parse_cli(int argc, char **argv, params_t &params){
    bool got_input = false;
    bool new_counts_name = false;
    bool new_fasta_name = false;
    int c;
    while((c = getopt(argc, argv, "i:k:vo:D:dhe:s:x:c:pP")) != -1){
        switch(c){
            case 'i':
                params.input_file_name = string(optarg);
                got_input = true;
                break;
            case 'o':
                params.fasta_file_name = string(optarg);
                new_fasta_name = true;
                break;
            case 'c':
                params.counts_file_name = string(optarg);
                new_counts_name = true;
                break;
            case 'k':
                params.kmer_size = atoi(optarg);
                if(params.kmer_size <= 0) {
                    cerr << "parse_cli(): Need a positive kmer size!" << endl;
                    exit(EXIT_FAILURE);
                }
                if(params.kmer_size % 2 == 0){
                    cerr << "parse_cli(): You should use an odd kmer size in order to avoid auto-loops in the DBG!" << endl;
                    exit(EXIT_SUCCESS);
                }
                break;
            case 'D':
                params.visited_depth = atoi(optarg);
                if(params.visited_depth <= 0) {
                    cerr << "parse_cli(): Need a positive depth!" << endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'p':
                params.phases = true;
                break;
            case 'P':
                params.profile = true;
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

    // check for input file
    if(!got_input){
        print_help(params);
        exit(EXIT_FAILURE);
    }

    // --- derive names ---
    // input = "../experiments/SRR001665_1.unitigs.fa"
    // get a base name removing BCALM extension ".unitigs.fa"
    auto ext_pos = params.input_file_name.rfind(".unitigs.fa");
    auto slash_pos = params.input_file_name.rfind('/');
    auto name_pos = (slash_pos == string::npos) ? 0 : slash_pos + 1; // if file is in PWD start from 0
    auto base_name = params.input_file_name.substr(name_pos, ext_pos - name_pos);

    if(!new_fasta_name)
        params.fasta_file_name = base_name + ".ustar.fa";
    if(!new_counts_name)
        params.counts_file_name = base_name + ".ustar" + encoding_suffixes.at(params.encoding) + ".counts";
}

int main(int argc, char **argv) {
    cout << "===== Unitig STitch Advanced constRuction (USTAR) v" << VERSION << " =====\n";
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
    if(params.debug) {
        if(!dbg.verify_input()){
            cerr << "Bad input file" << endl;
            exit(EXIT_FAILURE);
        }
    }

    if(params.profile)
        exit(0);

    // choose SPSS sorter
    Sorter sorter(params.seeding_method, params.extending_method, params.debug);
    // make an SPSS
    SPSS spss(&dbg, &sorter, params.visited_depth , params.phases, params.debug);

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
