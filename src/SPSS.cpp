//
// Created by enrico on 21/12/22.
//

#include <iostream>
#include "SPSS.h"

SPSS::SPSS(DBG *dbg){
    this->dbg = dbg;
}

void SPSS::simpler_test(const vector<size_t> &nodes, const vector<bool> &forwards) {
    if(dbg->check_path_consistency(nodes, forwards))
        cout << "Path is consistent!" << endl;
    else
        cout << "Path is NOT consistent!" << endl;

    string contig = dbg->spell(nodes, forwards);
    cout << contig << " (" << contig.length() << ")" << endl;

    vector<uint32_t> counts;
    dbg->get_counts(nodes, forwards, counts);
    for(auto count : counts)
        cout << count << " ";
    cout << endl;
}
