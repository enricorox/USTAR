//
// Created by enrico on 20/12/22.
//

#ifndef USTAR_DBG_H
#define USTAR_DBG_H

#include <string>
#include <vector>

#define MAX_LINE_LEN 100000

using namespace std;

struct arc_t{
    uint32_t successor;
    bool forward;
    bool to_forward;
};

struct node_t{
    uint32_t length;
    uint32_t median_abundance;
    double mean_abundance;
    string unitig;
    vector<uint32_t> abundances;
    vector<arc_t> arcs;
};

class DBG{
    string bcalm_file_name;
    uint32_t kmer_size = 0;
    vector<node_t> nodes;
    size_t n_arcs = 0;
    size_t n_kmers = 0;
    double avg_unitig_len = 0;
    double avg_abundances = 0;

    /**
     * Parse the BCALM2 file
     * @param bcalm_file_name Path of a file made by BCALM2 with '-all-abundance-counts'
     */
    void parse_bcalm_file();

    /**
     * Check wether two adjacent node labels overlap
     * @param node a dBG node
     * @param arcs a dBG arc
     * @return
     */
    bool overlaps(const node_t &node, const arc_t &arcs);

public:
    /**
     * Construct a de Bruijn Graph from a BCALM2 file
     * @param bcalm_file_name
     * @param kmer_size
     */
    DBG(const string &bcalm_file_name, uint32_t kmer_size);

    ~DBG();

    /**
     * Print some dBG stats
     */
    void print_info();

    /**
     * Verify that if there is an arcs between two nodes then they share a k-1 substring
     * @return true if all nodes satisfies that condition
     */
    bool verify_overlaps();

    /**
     * Write a BCALM2 like file from dBG in memory
     * @param file_name fasta file name
     */
    void to_bcalm_file(const string &file_name);

    /**
     * Compare BCALM2 file with to_bcalm_file()
     * @return true if the files are the same (spaces removed)
     */
    bool validate();

    /**
     * Compute the reverse complement
     * @param s a nucleotide sequence
     * @return the reverse-complement of s
     */
    static string reverse_complement(const string &s);

    /**
     * Get nodes reachable from node
     * @param node the current node ID
     * @param to_nodes a list of node IDs
     * @param to_forwards how a node must be read
     * @param mask nodes to filter
     */
    void get_nodes_from(uint32_t node, vector<bool> &forwards, vector<size_t> &to_nodes, vector<bool> &to_forwards, const vector<bool> &mask);

    /**
     * Like get_nodes_from() but nodes are path-consistent
     * @param node the current node ID
     * @param forward the direction of the current node
     * @param to_nodes a list of node IDs
     * @param to_forwards how a node must be read
     * @param mask nodes to filter
     */
    void get_consistent_nodes_from(uint32_t node, bool forward, vector<size_t> &to_nodes, vector<bool> &to_forwards, const vector<bool> &mask);

    /**
     * Compute the spell of this path, gluing node labels
     * @param path_nodes the path nodes
     * @param forwards how nodes must be read
     * @return the spell of the path
     */
    string spell(const vector<size_t> &path_nodes, const vector<bool> &forwards);

    /**
     * Check whether the path is path-consistent
     * @param path_nodes the nodes of the path
     * @param forwards how nodes must be read
     * @return true if the path is path-consistent
     */
    bool check_path_consistency(const vector<size_t> &path_nodes, const vector<bool> &forwards);

    /**
     * Get the kmers abundances in the given path
     * @param path_nodes the nodes of the path
     * @param forwards how nodes must be read
     * @param counts abundances are returned here
     */
    void get_counts(const vector<size_t> &path_nodes, const vector<bool> &forwards, vector<size_t> &counts);

    uint32_t get_n_kmers() const;

    uint32_t get_n_nodes() const;

    bool verify_input();

    const node_t &get_node(size_t node);

    uint32_t get_kmer_size();

    const vector<node_t> *get_nodes();
};

#endif //USTAR_DBG_H
