#include <iostream>
#include <getopt.h>
#include "DBG.h"
#include "SPSS.h"
#include "Encoder.h"

#define VERSION "0.1"

using namespace std;

struct params_t{
    string input_file{};
    string output_file = "out";
    string fasta_file_name = output_file + ".ustar.fa";
    string counts_file_name = output_file + ".ustar.counts";

    int kmer_size = 31;

    bool debug = false;

    encoding_t encoding = encoding_t::PLAIN;
    seeding_method_t seeding_method = seeding_method_t::FIRST;
    extending_method_t extending_method = extending_method_t::FIRST;
};

const map<encoding_t, string> encoding_suffixes = {
        {encoding_t::PLAIN, ""},
        {encoding_t::RLE, ".rle"},
        {encoding_t::FLIP, ".flip"},
        {encoding_t::AVG_RLE, ".avg_rle"},
        {encoding_t::FLIP_RLE, ".flip_rle"},
        {encoding_t::AVG_FLIP_RLE, ".avg_flip_rle"},
        {encoding_t::BINARY, ".bin"}
};

const map<string, encoding_t> encoding_names = {
        {"plain", encoding_t::PLAIN},
        {"flip", encoding_t::FLIP},
        {"rle", encoding_t::RLE},
        {"avg_rle", encoding_t::AVG_RLE},
        {"flip_rle", encoding_t::FLIP_RLE},
        {"avg_flip_rle", encoding_t::AVG_FLIP_RLE},
        {"bin", encoding_t::BINARY}
};

const map<string, seeding_method_t> seeding_method_names = {
        {"f", seeding_method_t::FIRST},
        {"-ma", seeding_method_t::LOWER_MEDIAN_ABUNDANCE},
        {"-aa", seeding_method_t::LOWER_AVERAGE_ABUNDANCE},
        {"=a", seeding_method_t::SIMILAR_ABUNDANCE},
        {"-l", seeding_method_t::SMALLER_LENGTH},
        {"+l", seeding_method_t::BIGGER_LENGTH},
        {"-c", seeding_method_t::LESS_CONNECTED},
        {"+c", seeding_method_t::MORE_CONNECTED}
};

const map<string, extending_method_t> extending_method_names = {
        {"f", extending_method_t::FIRST},
        {"=a", extending_method_t::SIMILAR_ABUNDANCE},
        {"-l", extending_method_t::SMALLER_LENGTH},
        {"+l", extending_method_t::BIGGER_LENGTH},
        {"-c", extending_method_t::LESS_CONNECTED},
        {"+c", extending_method_t::MORE_CONNECTED}
};

template<typename T>
string inv_map(const map<string, T> &m, const T &name){
    for(auto &p : m)
        if(p.second == name)
            return p.first;
    return "?";
}

void print_help(const params_t &params){
    cout << "Find a Spectrum Preserving String Set (aka simplitigs) for the input file.\n";
    cout << "Compute the kmer counts vector.\n\n";

    cout << "Usage: ./USTAR -i <bcalm_file>\n\n";
    cout << "Options:\n";

    cout << "   -k  kmer size, must be the same of BCALM2 [" << params.kmer_size << "]\n\n";

    cout << "   -o  output file base name [" << params.output_file << "]\n\n";

    cout << "   -s  seeding method [" << inv_map<seeding_method_t>(seeding_method_names, params.seeding_method) << "]\n";
    cout << "       f               choose the first seed available\n";
    cout << "       -ma             choose the seed with lower median abundance\n";
    cout << "       -aa             choose the seed with lower average abundance\n";
    cout << "       =a              choose the seed with most similar abundance to the last used node\n";
    cout << "       -l              choose the seed with smaller length\n";
    cout << "       +l              choose the seed with bigger length\n";
    cout << "\n";

    cout << "   -x  extending method [" << inv_map<extending_method_t>(extending_method_names, params.extending_method) << "]\n";
    cout << "       f               choose the first successor available\n";
    cout << "       =a              choose the successor with most similar abundance to the last used node\n";
    cout << "       -l              choose the successor with smaller length\n";
    cout << "       +l              choose the successor with bigger length\n";
    cout << "\n";

    cout << "   -e  encoding [" << inv_map<encoding_t>(encoding_names, params.encoding)<< "]\n";
    cout << "       plain           do not use any encoding\n";
    cout << "       rle             use special Run Length Encoding\n";
    cout << "       avg_rle         sort simplitigs by average counts and use RLE\n";
    cout << "       flip_rle        make contiguous runs by flipping simplitigs if necessary and use RLE\n";
    cout << "       avg_flip_rle    make contiguous runs by sorting by average, flipping simplitigs if necessary and use RLE\n";
    cout << "\n";

    cout << "   -d  debug [" << (params.debug?"true":"false") << "]\n\n";

    cout << "   -v  print version and author\n\n";

    cout << "   -h  print this help\n\n" << endl;
}

void print_params(const params_t &params){
    cout << "Params:\n";
    cout << "   input file:             " << params.input_file << "\n";
    cout << "   kmer size:              " << params.kmer_size << "\n";
    cout << "   output file base name:  " << params.output_file << "\n";
    cout << "   seeding method:         " << inv_map<seeding_method_t>(seeding_method_names, params.seeding_method) << "\n";
    cout << "   extending method:       " << inv_map<extending_method_t>(extending_method_names, params.extending_method) << "\n";
    cout << "   encoding:               " << inv_map<encoding_t>(encoding_names, params.encoding) << "\n";
    cout << "   debug:                  " << (params.debug?"true":"false") << "\n";
    cout << endl;
}

void parse_cli(int argc, char **argv, params_t &params){
    bool done = false;
    int c;
    while((c = getopt(argc, argv, "i:k:vo:dhe:s:x:")) != -1){
        switch(c){
            case 'i':
                if (!optarg) {
                    cerr << "parse_cli(): Need an input file name!" << endl;
                    exit(EXIT_FAILURE);
                }
                params.input_file = string(optarg);
                done = true;
                break;
            case 'o':
                if (!optarg) {
                    cerr << "parse_cli(): Need an output file name!" << endl;
                    exit(EXIT_FAILURE);
                }
                params.output_file = string(optarg);
                params.fasta_file_name = params.output_file + ".ustar.fa";
                params.counts_file_name =
                        params.output_file + ".ustar.counts" + encoding_suffixes.at(params.encoding);
                break;
            case 'k':
                if(optarg)
                    params.kmer_size = atoi(optarg);
                if(!optarg || params.kmer_size <= 0) {
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
                if (!optarg) {
                    cerr << "parse_cli(): need a method for encoding!" << endl;
                    exit(EXIT_FAILURE);
                }
                // is a valid encoding?
                if(encoding_names.find(optarg) == encoding_names.end()){
                    cerr << "parse_cli(): " << optarg << " is not a valid encoding" <<endl;
                    exit(EXIT_FAILURE);
                }
                params.encoding = encoding_names.at(optarg);
                params.counts_file_name = params.output_file
                                          + ".ustar.counts" + encoding_suffixes.at(params.encoding);
                break;
            case 's': // seed method
                if (!optarg) {
                    cerr << "parse_cli(): need a method for seeding" << endl;
                    exit(EXIT_FAILURE);
                }
                if(seeding_method_names.find(optarg) == seeding_method_names.end()){
                    cerr << "parse_cli(): " << optarg << " is not a valid seed method" <<endl;
                    exit(EXIT_FAILURE);
                }
                params.seeding_method = seeding_method_names.at(optarg);
                break;
            case 'x': // extension method
                if (!optarg) {
                    cerr << "parse_cli(): need a method for extension" << endl;
                    exit(EXIT_FAILURE);
                }
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
                cerr << "parse_cli(): missing argument\n\n";
                print_help(params);
                exit(EXIT_FAILURE);
            default:
                cerr << "parse_cli(): unknown parameter in optstring '" << c << "'\n\n";
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
    cout << "===== Unitig STitch STar (USTAR) v" VERSION " =====" << endl;
    // cli parameters
    params_t params;
    parse_cli(argc, argv, params);
    print_params(params);

    // make a dBG
    cout << "Reading the input file..." << endl;
    DBG dbg(params.input_file, params.kmer_size);
    dbg.print_stat();

    // verify input
    if(params.debug)
        dbg.verify_input();

    // choose SPSS sorter
    Sorter sorter(params.seeding_method, params.extending_method, params.debug);
    // make an SPSS
    SPSS spss(&dbg, &sorter, params.debug);

    // compute simplitigs
    cout << "Computing a path cover..." << endl;
    spss.compute_path_cover();
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
