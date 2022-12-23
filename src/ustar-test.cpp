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
            {"../unit-test/test1",
             "../unit-test/test2",
             "../unit-test/test3",
             "../unit-test/test4",
             "../unit-test/test5",
             "../unit-test/test6",
             "../unit-test/test7",
             "../unit-test/test8",
             "../unit-test/test9",
             "../unit-test/test10",
             "../unit-test/test11"
            };

    for(auto &file_name_base : file_names){
        string file_name = file_name_base + ".unitigs.fa";
        cout << "Test file: " << file_name << endl;

        DBG dbg(file_name, 3);
        dbg.verify_input();
        dbg.print_info();

        Sorter sorter;
        SPSS spss(&dbg, &sorter, true);
        spss.compute_path_cover();
        spss.extract_simplitigs_and_counts();
        spss.print_info();
        Encoder encoder(spss.get_simplitigs(), spss.get_counts());
        encoder.to_fasta_file(file_name_base + ".ustar.fa");
        encoder.to_counts_file(file_name_base + ".ustar.counts");

        // verify output
        ifstream fasta;
        fasta.open(file_name_base + ".ustar.fa");
        string result, correct_result = "ACTGG";
        fasta >> result >> result; // ignore ">"
        cout << "Result: " << result << endl;
        if(result != correct_result && result != DBG::reverse_complement(correct_result)){
            cerr << "Output must be: ACTGG or " << DBG::reverse_complement(correct_result) << "!" << endl;
            cerr << "Found " << result << " instead." << endl;
            exit(EXIT_FAILURE);
        }else
            cout << "YES! Output is correct!" << endl;
        cout << endl;
        fasta.close();
    }
}