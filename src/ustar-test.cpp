//
// Created by enrico on 22/12/22.
//
#include <iostream>
#include "SPSS.h"

using namespace std;

int main(){
    cout << "===== USTAR unit test =====" << endl;
    vector<string> file_names =
            {"../unit-test/test1",
             "../unit-test/test2",
             "../unit-test/test3",
             "../unit-test/test4"
            };

    for(auto &file_name_base : file_names){
        string file_name = file_name_base + ".unitigs.fa";
        cout << "Test file: " << file_name << endl;
        DBG dbg(file_name, 3);
        dbg.verify_input();
        dbg.print_info();
        SPSS spss(&dbg);
        spss.extract_simplitigs();
        spss.print_info();
        spss.to_fasta_file(file_name_base + ".ustar.fa");
        cout << endl;
    }
}