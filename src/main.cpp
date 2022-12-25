#include <iostream>
#include <getopt.h>
#include "DBG.h"
#include "SPSS.h"
#include "Encoder.h"

#define VERSION "0.1"

using namespace std;

struct params_t{
    string input_file = "../experiments/k31.a1.unitigs.fa";
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
        {"=a", seeding_method_t::SIMILAR_ABUNDANCE}
};

const map<string, extending_method_t> extending_method_names = {
        {"f", extending_method_t::FIRST},
        {"=a", extending_method_t::SIMILAR_ABUNDANCE}
};

template <typename T>
string inv_map(const map<string, T> m, T name){
    for(auto &p : m)
        if(p.second == name)
            return p.first;
    return "?";
}

void print_help(const params_t &params){
    cout << "Find a Spectrum Preserving String Set (aka simplitigs) for the input file.\n";
    cout << "Compute the kmer simplitigs_counts vector.\n\n";

    cout << "Usage: ./USTAR -i <bcalm_file>\n\n";
    cout << "Options:\n";

    cout << "\t-k \tkmer size [" << params.kmer_size << "]\n\n";

    cout << "\t-o \toutput file base name [" << params.output_file << "]\n\n";

    cout << "\t-s \tseeding method [" << inv_map<seeding_method_t>(seeding_method_names, params.seeding_method) << "]\n";
    cout << "\t\tf\tchoose the first seed available\n";
    cout << "\t\t-ma\tchoose the seed with lower median abundance\n";
    cout << "\t\t=a\tchoose the seed with most similar abundance to the last used node\t";
    cout << "\n";

    cout << "\t-x \textending method [" << inv_map<extending_method_t>(extending_method_names, params.extending_method) << "]\n";
    cout << "\t\tf\tchoose the first successor available\n";
    cout << "\t\t=a\tchoose the successor with most similar abundance to the last used node\t";
    cout << "\n";

    cout << "\t-e \tencoding [" << inv_map<encoding_t>(encoding_names, params.encoding)<< "]\n";
    cout << "\t\tplain\tdo not use any encoding\n";
    cout << "\t\trle\tuse special Run Length Encoding\n";
    cout << "\t\tavg_rle\tsort simplitigs by average counts and use RLE\n";
    cout << "\t\tflip_rle\tmake contiguous runs by flipping simplitigs if necessary and use RLE\n";
    cout << "\n";

    cout << "\t-d \tdebug [" << (params.debug?"true":"false") << "]\n\n";

    cout << "\t-v \tprint version and author\n\n";

    cout << "\t-h \tprint this help\n\n" << endl;
}

void print_params(const params_t &params){
    cout << "Params:\n";
    cout << "\tinput file: " << params.input_file << "\n";
    cout << "\tkmer size: " << params.kmer_size << "\n";
    cout << "\toutput file base name: " << params.output_file << "\n";
    cout << "\tseeding method: " << inv_map<seeding_method_t>(seeding_method_names, params.seeding_method) << "\n";
    cout << "\textending method: " << inv_map<extending_method_t>(extending_method_names, params.extending_method) << "\n";
    cout << "\tencoding: " << inv_map<encoding_t>(encoding_names, params.encoding) << "\n";
    cout << "\tdebug: " << (params.debug?"true":"false") << "\n";
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
                // if(encoding_names.find(optarg) == encoding_names.end()){}
                params.encoding = encoding_names.at(optarg);
                params.counts_file_name = params.output_file
                                          + ".ustar.counts" + encoding_suffixes.at(params.encoding);
                break;
            case 's': // seed method
                if (!optarg) {
                    cerr << "parse_cli(): need a method for seeding" << endl;
                    exit(EXIT_FAILURE);
                }
                params.seeding_method = seeding_method_names.at(optarg);
                break;
            case 'x': // extension method
                if (!optarg) {
                    cerr << "parse_cli(): need a method for extension" << endl;
                    exit(EXIT_FAILURE);
                }
                params.extending_method = extending_method_names.at(optarg);
                break;
            case 'h':
                print_help(params);
                exit(EXIT_SUCCESS);
            default:
                cerr << "parse_cli(): unknown parameter -" << c << endl;
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
    cout << "Extracting simplitigs and kmers simplitigs_counts..." << endl;
    spss.extract_simplitigs_and_counts();
    spss.print_stat();

    Encoder encoder(spss.get_simplitigs(), spss.get_counts(), params.debug);
    encoder.encode(params.encoding);
    encoder.print_stat();
    encoder.to_fasta_file(params.fasta_file_name);
    cout << "Simplitigs written to disk: " << params.fasta_file_name << endl;
    encoder.to_counts_file(params.counts_file_name);
    cout << "Counts written to disk: " << params.counts_file_name << endl;

    return EXIT_SUCCESS;
}
