//
// Created by enrico on 22/12/22.
//
#include <iostream>
#include <fstream>
#include "SPSS.h"

using namespace std;

int main(){
    cout << "===== USTAR unit test =====" << endl;
    vector<string> file_names =
            {"../tests/test1",
             "../tests/test2",
             "../tests/test3",
             "../tests/test4",
             "../tests/test5",
             "../tests/test6",
             "../tests/test7",
             "../tests/test8",
             "../tests/test9",
             "../tests/test10",
             "../tests/test11"
            };

    for(auto &file_name_base : file_names){
        string file_name = file_name_base + ".unitigs.fa";
        cout << "Test file: " << file_name << endl;

        DBG dbg(file_name, 3);
        dbg.verify_input();
        dbg.print_stat();

        Sorter sorter(seeding_method_t::HIGHER_AVERAGE_ABUNDANCE, extending_method_t::LESS_CONNECTED);
        SPSS spss(&dbg, &sorter, false, true);
        spss.compute_path_cover();
        spss.extract_simplitigs_and_counts();
        spss.print_stats();
        Encoder encoder(spss.get_simplitigs(), spss.get_counts());
        encoder.to_fasta_file(file_name_base + ".ustar.fa");
        encoder.to_counts_file(file_name_base + ".ustar.simplitigs_counts");

        // verify output
        ifstream fasta;
        fasta.open(file_name_base + ".ustar.fa");
        string result, correct_result = "ACTGG";
        fasta >> result >> result; // ignore ">"
        cout << "Result: " << result << endl;
        if(result != correct_result && result != DBG::reverse_complement(correct_result)){
            cerr << "Output must be: "<< correct_result <<" or " << DBG::reverse_complement(correct_result) << "!" << endl;
            cerr << "Found " << result << " instead." << endl;
            exit(EXIT_FAILURE);
        }else
            cout << "YES! Output is correct!" << endl;
        cout << endl;
        fasta.close();
    }
}