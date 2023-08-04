//
// Created by enrico on 21/12/22.
//

#include <iostream>
#include <stack>
#include <fstream>
#include <algorithm>
#include <set>
#include "SPSS.h"

SPSS::SPSS(DBG *dbg, Sorter *sorter, int depth, bool phases, bool debug){
    this->debug = debug;
    this->depth = depth;
    this->phases = phases;

    this->dbg = dbg;
    n_nodes = dbg->get_n_nodes();
    visited.resize(n_nodes, false);

    this->sorter = sorter;
    sorter->init(dbg->get_nodes(), &visited);
}

const vector<string> * SPSS::get_simplitigs(){
    if(simplitigs.empty()){
        cerr << "There are no simplitigs!" << endl;
        exit(EXIT_FAILURE);
    }
    return &simplitigs;
}

const vector<vector<uint32_t>> * SPSS::get_counts(){
    if(simplitigs.empty()){
        cerr << "There are no simplitigs!" << endl;
        exit(EXIT_FAILURE);
    }
    return &counts;
}

size_t SPSS::compute_contig_length(const vector<node_idx_t>& nodes){
    size_t len = 0;
    for(auto n: nodes){
        len += dbg->get_node(n).length;
    }
    len = len - (dbg->get_kmer_size() - 1) * (nodes.size() - 1);
    if(debug) cout << "contig length: " << len << "\n";
    return len;
}

void SPSS::jump_visited(node_idx_t node, bool forward, size_t max_depth, vector<node_idx_t> &nodes, vector<bool> &forwards){
    // --- is there a node with non-visited neighbours? ---
    // --- find the one with minimum connectivity ---
    vector<bool> dummy_visited;
    vector<node_idx_t> curr_to_nodes; vector<bool> curr_to_forwards;

    // visited node in BFS
    set<node_idx_t> visited_bfs;
    visited_bfs.insert(node);

    // paths to each node
    map<node_idx_t, vector<node_idx_t>> paths_nodes;
    map<node_idx_t, vector<bool>> paths_forwards;

    // path to starting node
    paths_nodes[node] = vector<node_idx_t>();
    paths_nodes[node].push_back(node);
    paths_forwards[node] = vector<bool>();
    paths_forwards[node].push_back(forward);

    deque<node_idx_t> q_nodes;
    deque<bool> q_forwards;

    q_nodes.push_back(node);
    q_forwards.push_back(forward);

    size_t z_min = UINT64_MAX;
    node_idx_t node_min = node;
    while(!q_nodes.empty()){
        // pop node & orientation
        node_idx_t curr = q_nodes.front();
        q_nodes.pop_front();
        bool curr_forward = q_forwards.front();
        q_forwards.pop_front();

        // contig length check
        if(paths_nodes.at(curr).size() > 2) {
            auto doubled = vector<node_idx_t>(paths_nodes.at(curr).begin() + 1,
                                              paths_nodes.at(curr).end() - 1);
            if (compute_contig_length(doubled) > 2 * dbg->get_kmer_size() - 2) continue;
        }

        // push its children to q and save path
        dbg->get_consistent_nodes_from(curr, curr_forward, curr_to_nodes, curr_to_forwards, dummy_visited);
        for(size_t i = 0; i < curr_to_nodes.size(); i++){
            // find loops
            //if(paths_nodes.find(curr_to_nodes[i]) != paths_nodes.end()) continue;
            if(visited_bfs.find(curr_to_nodes[i]) != visited_bfs.end()) continue;
            visited_bfs.insert(curr_to_nodes[i]);

            q_nodes.push_back(curr_to_nodes[i]);
            q_forwards.push_back(curr_to_forwards[i]);

            // save path nodes
            vector<node_idx_t> new_path_nodes = vector<node_idx_t>(paths_nodes.at(curr));
            new_path_nodes.push_back(curr_to_nodes[i]);
            paths_nodes[curr_to_nodes[i]] =  new_path_nodes;
            // save path forwards
            vector<bool> new_path_forwards = vector<bool>(paths_forwards.at(curr));
            new_path_forwards.push_back(curr_to_forwards[i]);
            paths_forwards[curr_to_nodes[i]] = new_path_forwards;
        }

        // check depth
        if(debug) cout << "size(" << curr << ") = " << paths_nodes.at(curr).size() << "\n";
        if(paths_nodes.at(curr).size() > max_depth + 2) {
            if(debug) cout << "too deep: break\n";
            break;
        }

        // check unvisited
        if(!visited.at(curr)) {
            max_depth = paths_nodes.at(curr).size() - 2;

            // accumulate minimum
            if (dbg->get_node(curr).arcs.size() < z_min) {
            //if (dbg->get_node(curr).length < z_min) {
                z_min = dbg->get_node(curr).arcs.size();
                //z_min = dbg->get_node(curr).length;
                node_min = curr;
            }
        }
    }
    // return minimum cost path
    nodes = paths_nodes.at(node_min);
    forwards = paths_forwards.at(node_min);
}

size_t SPSS::jump_visited_phase(node_idx_t node, bool forward, size_t max_depth, vector<node_idx_t> &nodes, vector<bool> &forwards){
    // --- is there a node with non-visited neighbours? ---
    // --- find the one with minimum connectivity ---
    vector<bool> dummy_visited;
    vector<node_idx_t> child_to_nodes; vector<bool> child_to_forwards;

    // paths to each node
    map<node_idx_t, vector<node_idx_t>> paths_nodes;
    map<node_idx_t, vector<bool>> paths_forwards;

    // path to starting node
    paths_nodes[node] = vector<node_idx_t>();
    paths_nodes[node].push_back(node);
    paths_forwards[node] = vector<bool>();
    paths_forwards[node].push_back(forward);

    deque<node_idx_t> q_nodes;
    deque<bool> q_forwards;

    q_nodes.push_back(node);
    q_forwards.push_back(forward);

    size_t z_min = UINT64_MAX;
    node_idx_t node_min = node;
    while(!q_nodes.empty()){
        // pop node & orientation
        node_idx_t child = q_nodes.front();
        q_nodes.pop_front();
        bool child_forward = q_forwards.front();
        q_forwards.pop_front();

        // push its children to q and save path
        dbg->get_consistent_nodes_from(child, child_forward, child_to_nodes, child_to_forwards, dummy_visited);
        for(size_t i = 0; i < child_to_nodes.size(); i++){
            // find and break loops
            if(paths_nodes.find(child_to_nodes[i]) != paths_nodes.end()) continue;

            q_nodes.push_back(child_to_nodes[i]);
            q_forwards.push_back(child_to_forwards[i]);

            // save path nodes
            vector<node_idx_t> new_path_nodes = vector<node_idx_t>(paths_nodes.at(child));
            new_path_nodes.push_back(child_to_nodes[i]);
            paths_nodes[child_to_nodes[i]] =  new_path_nodes;
            // save path forwards
            vector<bool> new_path_forwards = vector<bool>(paths_forwards.at(child));
            new_path_forwards.push_back(child_to_forwards[i]);
            paths_forwards[child_to_nodes[i]] = new_path_forwards;
        }

        // check depth
        if(debug) cout << "size(" << child << ") = " << paths_nodes.at(child).size() << "\n";
        if(paths_nodes.at(child).size() - 1 > max_depth) {
            if(debug) cout << "too deep: break\n";
            break;
        }

        // check if we can link paths
        if(heads.find(child) != heads.end()){
            max_depth = paths_nodes.at(child).size() - 1;

            // accumulate minimum
            if (dbg->get_node(child).arcs.size() < z_min) {
                z_min = dbg->get_node(child).arcs.size();
                node_min = child;
            }
        }
    }
    // return minimum cost path
    nodes = paths_nodes.at(node_min);
    forwards = paths_forwards.at(node_min);
    return node_min;
}

void SPSS::extends(vector<node_idx_t> &path_nodes, vector<bool> &path_forwards) {
    vector<node_idx_t> to_nodes; vector<bool> to_forwards; vector<bool> forwards;

    // add the seed
    node_idx_t seed = path_nodes.back();
    bool forward = path_forwards.back();
    visited.at(seed) = true;

    // extension vars
    bool to_forward;
    size_t node, successor;

    // ----- extending -----
    node = seed;
    while (true) {
        // explore
        dbg->get_consistent_nodes_from(node, forward, to_nodes, to_forwards, visited);

        if (to_nodes.empty()) {
            if(phases || (depth <= 0)) break;

            if(debug) cout << "No more unvisited; jumping...\n";
            vector<node_idx_t> jump_nodes; vector<bool> jump_forwards;
            jump_visited(node, forward, depth, jump_nodes, jump_forwards);
            if(jump_nodes.size() < 2) {
                if(debug) cout << "No node found; break and find a new seed\n";
                break;
            }
            path_nodes.insert(path_nodes.end(), jump_nodes.begin() + 1, jump_nodes.end());
            path_forwards.insert(path_forwards.end(), jump_forwards.begin() + 1, jump_forwards.end());
            node = jump_nodes.back();
            forward = jump_forwards.back();
            visited.at(node) = true;
            continue;
        }
        successor = sorter->next_successor(node, forward, to_nodes, to_forwards, to_forward);

        // update
        node = successor;
        forward = to_forward;
        visited.at(node) = true;
        path_nodes.push_back(node);
        path_forwards.push_back(forward);
    }

    if(debug && !dbg->check_path_consistency(path_nodes, path_forwards))
        cerr << "Inconsistent path!\n";
}

void reverse_path(vector<node_idx_t> &path_nodes, vector<bool> &path_forwards){
    std::reverse(path_nodes.begin(), path_nodes.end());
    std::reverse(path_forwards.begin(), path_forwards.end());
    for(auto && forward : path_forwards)
        forward = (!forward);
}

void SPSS::compute_path_cover(bool two_way) {
    // reset path cover if already computed
    path_cover_nodes.clear();
    path_cover_forwards.clear();
    std::fill(visited.begin(), visited.end(), false);

    vector<node_idx_t> path_nodes; vector<bool> path_forwards;
    while(sorter->has_seed()){
        path_nodes.clear(); path_forwards.clear();

        path_nodes.push_back(sorter->next_seed());
        path_forwards.push_back(true);

        // forward extending
        extends(path_nodes, path_forwards);
        if(two_way) {
            if(debug) cout << "Doing two-way...\n";
            // backward extending
            reverse_path(path_nodes, path_forwards);
            extends(path_nodes, path_forwards);
            if(debug) cout << "Done.\n";
        }

        path_cover_nodes.push_back(path_nodes);
        path_cover_forwards.push_back(path_forwards);

        if(phases) {
            // put starts & ends in hash map
            heads[path_nodes.front()] = path_nodes.size() - 1;
            tails[path_nodes.back()] = path_nodes.size() - 1;
        }
    }
    if(phases) {
        // --- find link between path using jumps ---
        // ...
    }
}

void SPSS::extract_simplitigs_and_counts(){
    n_simplitigs = path_cover_nodes.size();
    if(n_simplitigs == 0){
        cerr << "extract_simplitigs_and_counts(): Need to compute a path cover first!" << endl;
        exit(EXIT_FAILURE);
    }

    vector<uint32_t> simplitig_counts;
    counts.reserve(n_simplitigs);
    for(size_t i = 0; i < n_simplitigs; i++){
        // extract simplitigs
        simplitigs.push_back(dbg->spell(path_cover_nodes[i], path_cover_forwards[i]));

        // extract simplitigs_counts
        simplitig_counts.clear();
        dbg->get_counts(path_cover_nodes[i], path_cover_forwards[i], simplitig_counts);
        counts.push_back(simplitig_counts);
    }
}

void SPSS::print_stats(){
    if(n_simplitigs == 0){
        cerr << "print_stats(): Need to extract simplitigs first!" << endl;
        exit(EXIT_FAILURE);
    }
    size_t c_length = 0;
    for(auto &path : path_cover_nodes){
        c_length += dbg->get_node(path[0]).length;
        for(size_t i = 1; i < path.size(); i++)
            c_length += dbg->get_node(path[i]).length - (dbg->get_kmer_size() - 1);
    }

    cout << "\n";
    cout << "Simplitigs stats:\n";
    cout << "   number of simplitigs (NS):              " << n_simplitigs << "\n";
    cout << "   cumulative length (CL):                 " << c_length << "\n";
    cout << "   nucleotide per kmer:                    " << (double) c_length / (double) dbg->get_n_kmers() << "\n";
    cout << "   average simplitigs length:              " << (double) c_length / (double) n_simplitigs << "\n";
    cout << "\n";
}

size_t SPSS::get_score() {
    return path_cover_nodes.size();
}
