//
// Created by enrico on 21/12/22.
//

#include <iostream>
#include "SPSS.h"

SPSS::SPSS(DBG *dbg){
    this->dbg = dbg;
    n_nodes = dbg->get_n_nodes();
    saturated.resize(n_nodes);
    fill(saturated.begin(), saturated.end(), false);
}

void SPSS::test() {
    vector<vector<size_t>> nodes_v = {{7, 157883, 101636}, {7, 157883, 82391}, {7, 157883, 114132}};
    vector<vector<bool>> forwards_v = {{true, true, true}, {true, true, true}, {true, true, true}};

    for(size_t i = 0; i < nodes_v.size(); i++) {
        auto &nodes = nodes_v.at(i);
        auto &forwards = forwards_v.at(i);

        if (dbg->check_path_consistency(nodes, forwards))
            cout << "Path is consistent!" << endl;
        else
            cout << "Path is NOT consistent!" << endl;

        string contig = dbg->spell(nodes, forwards);
        cout << contig << " (" << contig.length() << ")" << endl;

        vector<uint32_t> counts;
        dbg->get_counts(nodes, forwards, counts);
        for (auto count: counts)
            cout << count << " ";
        cout << endl;
    }
}
