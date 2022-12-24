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
    seeding_method_t seeding_method = seeding_method_t::DEFAULT_SEED;
    extending_method_t extending_method = extending_method_t::DEFAULT_EXTEND;
};

const map<encoding_t, string> encoding_suffixes = {
        {encoding_t::PLAIN, ""},
        {encoding_t::RLE, ".rle"},
        {encoding_t::AVG_RLE, ".avg_rle"}
};

const map<string, encoding_t> encoding_names = {
        {"d", encoding_t::PLAIN}, {"plain", encoding_t::PLAIN},
        {"rle", encoding_t::RLE},
        {"avg_rle", encoding_t::AVG_RLE}
};

const map<string, seeding_method_t> seeding_method_names = {
        {"d", seeding_method_t::DEFAULT_SEED},
};

const map<string, extending_method_t> extending_method_names = {
        {"d", extending_method_t::DEFAULT_EXTEND}
};

void print_help(){
    cout << "Find a Spectrum Preserving String Set (aka simplitigs) for the input file.\n";
    cout << "Compute the kmer simplitigs_counts vector.\n\n";
    cout << "Usage: ./USTAR -i <bcalm_file> -k <kmer_size>\n";
    cout << "Options:\n";
    cout << "\t-o \toutput file base name\n";
    cout << "\t-s \tseeding method\n";
    cout << "\t-x \textending method\n";
    cout << "\t-e \tencoding\n";
    cout << "\t-d \tdebug\n";
    cout << "\t-v \tprint version\n";
    cout << "\t-h \tprint this help\n" << endl;
}

void print_params(const params_t &params){
    cout << "Params:\n";
    cout << "\tinput file: " << params.input_file << "\n";
    cout << "\tkmer size: " << params.kmer_size << "\n";
    cout << "\toutput file base name: " << params.output_file << "\n";
    cout << "\tseeding method: "
         << [](seeding_method_t m){
                for(auto &p : seeding_method_names)
                    if(p.second == m)
                        return p.first;
                return string("?");
            }(params.seeding_method)
         << "\n";
    cout << "\textending method: "
         << [](extending_method_t m){
             for(auto &p : extending_method_names)
                 if(p.second == m)
                     return p.first;
             return string("?");
            }(params.extending_method)
         << "\n";
    cout << "\tencoding: "
         << [](encoding_t m){
             for(auto &p : encoding_names)
                 if(p.second == m)
                     return p.first;
             return string("?");
            }(params.encoding)
         << "\n";
    cout << "\tdebug: " << (params.debug?"true":"false") << "\n";
    cout << endl;
}

void parse_cli(int argc, char **argv, params_t &params){
    bool done = false;
    int c;
    while((c = getopt(argc, argv, "i:k:vo:dhe:s:x:")) != -1){
        done = true;
        switch(c){
            case 'i':
                if(optarg)
                    params.input_file = string(optarg);
                else {
                    cerr << "parse_cli(): Need an input file name!" << endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'o':
                if(optarg) {
                    params.output_file = string(optarg);
                    params.fasta_file_name = params.output_file + ".ustar.fa";
                    params.counts_file_name = params.output_file + ".ustar.count" + encoding_suffixes.at(params.encoding);
                }
                else {
                    cerr << "parse_cli(): Need an output file name!" << endl;
                    exit(EXIT_FAILURE);
                }
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
                cout << "Version " << VERSION << endl;
                exit(EXIT_SUCCESS);
            case 'd':
                params.debug = true;
                break;
            case 'e': // encoding
                if(optarg) {
                    params.encoding = encoding_names.at(optarg);
                    params.counts_file_name = params.output_file
                            + ".counts" + encoding_suffixes.at(params.encoding);
                }
                else{
                    cerr << "parse_cli(): need a method for encoding!" << endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 's': // seed method
                if (optarg)
                    params.seeding_method = seeding_method_names.at(optarg);
                else {
                    cerr << "parse_cli(): need a method for seeding" << endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'x': // extension method
                if (optarg)
                    params.extending_method = extending_method_names.at(optarg);
                else {
                    cerr << "parse_cli(): need a method for extension" << endl;
                    exit(EXIT_FAILURE);
                }
                break;
            case 'h':
                print_help();
                exit(EXIT_SUCCESS);
            default:
                cerr << "parse_cli(): unknown parameter -" << c << endl;
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
    Sorter sorter(params.seeding_method, params.extending_method);
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
