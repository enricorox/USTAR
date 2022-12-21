//
// Created by enrico on 20/12/22.
//

#ifndef USTAR_DBG_H
#define USTAR_DBG_H

#include <string>
#include <vector>

#define MAX_LINE_LEN 100000

using namespace std;

struct arcs_t{
    uint32_t successor;
    bool forward;
    bool to_forward;
};

struct node_t{
    uint32_t length;
    string unitig;
    double avg_abundance;
    vector<uint32_t> abundances;
    vector<arcs_t> arcs;
};

class DBG{
    string bcalm_file_name;
    int kmer_size = 0;
    vector<node_t> nodes;
    size_t n_arcs = 0;
    size_t n_kmers = 0;
    double avg_unitig_len = 0;
    double avg_abundances = 0;

    /**
     * Parse the BCALM2 file
     * @param bcalm_file_name
     */
    void parse_bcalm_file();

    bool overlaps(const node_t &node, const arcs_t &edge);

public:
    /**
     * Construct a de Bruijn Graph from a BCALM2 file
     * @param bcalm_file_name
     */
    DBG(const string &bcalm_file_name, int kmer_size);

    ~DBG();

    void print_info();

    bool verify_overlaps();

    void to_bcalm_file(const string &file_name);

    bool validate();

    static string reverse_complement(const string &s);

    void get_nodes_from(int node, vector<size_t> &to_nodes, vector<bool> &forwards);

    void get_consistent_nodes_from(int node, bool forward, vector<size_t> &to_nodes, vector<bool> &forwards);

    string spell(const vector<size_t> &path_nodes, const vector<bool> &forwards);

    bool check_path_consistency(const vector<size_t> &path_nodes, const vector<bool> &forwards);

    void get_counts(const vector<size_t> &path_nodes, const vector<bool> &forwards, vector<uint32_t> &counts);
};

#endif //USTAR_DBG_H
