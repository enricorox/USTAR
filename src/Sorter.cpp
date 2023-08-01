//
// Created by enrico on 22/12/22.
//

#include <iostream>
#include <algorithm>
#include <random>
#include "Sorter.h"
#include "commons.h"

const double EPSILON = 0.5;

Sorter::Sorter(seeding_method_t sorting_methods, extending_method_t extending_method, bool debug) { // NOLINT(cert-msc51-cpp)
    this->seeding_method = sorting_methods;
    this->extending_method = extending_method;
    this->debug = debug;
}

void Sorter::init(const vector<node_t> *dbg_nodes, const vector<bool> *spss_visited){
    this->visited = spss_visited;
    this->nodes = dbg_nodes;

    random_device rd;
    auto seed = rd();
    random_generator.seed(seed);
    if(debug)
        cout << "Random seed: " << seed << "\n";

    seed_order.reserve(dbg_nodes->size());
    for(size_t i = 0; i < dbg_nodes->size(); i++)
        seed_order.push_back(i);

    seed_index = 0;

    // sort!
    switch(seeding_method){
        case seeding_method_t::MORE_CONNECTED:{
            auto lambda = [this](size_t a, size_t b){
                return (*nodes)[a].arcs.size() > (*nodes)[b].arcs.size();
            };
            sort(seed_order.begin(), seed_order.end(), lambda);
            }
            break;
        case seeding_method_t::LESS_CONNECTED:{
            auto lambda = [this](size_t a, size_t b){
                return (*nodes)[a].arcs.size() < (*nodes)[b].arcs.size();
            };
            sort(seed_order.begin(), seed_order.end(), lambda);
            }
            break;
        case seeding_method_t::BIGGER_LENGTH: {
            auto lambda = [this](size_t a, size_t b){
                return nodes->at(a).length > nodes->at(b).length;
            };
            sort(seed_order.begin(), seed_order.end(), lambda);
            }
            break;
        case seeding_method_t::SMALLER_LENGTH:{
            auto lambda = [this](size_t a, size_t b){
                return nodes->at(a).length < nodes->at(b).length;
            };
            sort(seed_order.begin(), seed_order.end(), lambda);
            }
            break;
        case seeding_method_t::LOWER_MEDIAN_ABUNDANCE: {
            auto lambda = [this](size_t a, size_t b) {
                return nodes->at(a).median_abundance < nodes->at(b).median_abundance;
                //return nodes->at(a).abundances.front() + nodes->at(a).abundances.back()  < nodes->at(b).abundances.front() + nodes->at(b).abundances.back();
            };
            sort(seed_order.begin(), seed_order.end(), lambda);
            }
            break;
        case seeding_method_t::SIMILAR_ABUNDANCE:
            // no break here
        case seeding_method_t::LOWER_AVERAGE_ABUNDANCE: {
                auto lambda = [this](size_t a, size_t b) {
                    return nodes->at(a).average_abundance < nodes->at(b).average_abundance;
                };
                sort(seed_order.begin(), seed_order.end(), lambda);
            }
            break;
        case seeding_method_t::HIGHER_AVERAGE_ABUNDANCE: {
                auto lambda = [this](size_t a, size_t b) {
                    return nodes->at(a).average_abundance > nodes->at(b).average_abundance;
                };
                sort(seed_order.begin(), seed_order.end(), lambda);
            }
            break;
        case seeding_method_t::LESS_UNBALANCED: {
            auto lambda = [this](size_t a, size_t b) {
                uint32_t num_forward_a = 0; uint32_t num_backward_a = 0;
                uint32_t num_forward_b = 0; uint32_t num_backward_b = 0;
                for(auto arc : nodes->at(a).arcs)
                    if(arc.forward)
                        num_forward_a++;
                    else
                        num_backward_a++;
                for(auto arc : nodes->at(b).arcs)
                    if(arc.forward)
                        num_forward_b++;
                    else
                        num_backward_b++;
                return d(num_forward_a, num_backward_a) < d(num_forward_b, num_backward_b);
            };
            sort(seed_order.begin(), seed_order.end(), lambda);
        }
            break;
        case seeding_method_t::MORE_UNBALANCED: {
            auto lambda = [this](size_t a, size_t b) {
                uint32_t num_forward_a = 0; uint32_t num_backward_a = 0;
                uint32_t num_forward_b = 0; uint32_t num_backward_b = 0;
                for(auto arc : nodes->at(a).arcs)
                    if(arc.forward)
                        num_forward_a++;
                    else
                        num_backward_a++;
                for(auto arc : nodes->at(b).arcs)
                    if(arc.forward)
                        num_forward_b++;
                    else
                        num_backward_b++;
                return d(num_forward_a, num_backward_a) > d(num_forward_b, num_backward_b);
            };
            sort(seed_order.begin(), seed_order.end(), lambda);
        }
            break;
        case seeding_method_t::FIRST:
            break;
        case seeding_method_t::RANDOM:
            shuffle(seed_order.begin(), seed_order.end(), random_generator);
            break;
        default:
            cerr << "init(): unknown seeding method!" << endl;
            exit(EXIT_FAILURE);
    }
}

size_t Sorter::next_seed() {
    if(!has_seed()){
        cerr << "next_seed(): No seed available!" << endl;
        exit(EXIT_FAILURE);
    }
    if(seeding_method == seeding_method_t::SIMILAR_ABUNDANCE){
        if(first_node){
            first_node = false;
            return seed_order[seed_index];
        }
        size_t best = seed_index;
        auto d_best = d((*nodes)[last_node].median_abundance , (*nodes)[seed_order[best]].median_abundance);
        for(size_t i = seed_index; i < seed_order.size(); i++)
            if(!(*visited)[seed_order[i]]){
                auto d_i = d((*nodes)[last_node].median_abundance , (*nodes)[seed_order[i]].median_abundance);
                if(d_i < d_best)
                    best = i;
                if(d_i < EPSILON)
                    break;
            }
        swap(seed_order[seed_index], seed_order[best]);
    }
    last_node = seed_order[seed_index];
    return last_node;
}

bool Sorter::has_seed() {
    for(; seed_index < seed_order.size(); seed_index++)
        if(!(*visited)[seed_order[seed_index]])
            break;
    return seed_index < seed_order.size();
}

size_t Sorter::next_successor(node_idx_t seed, bool forward, vector<node_idx_t> &to_nodes, vector<bool> &to_forwards, bool &to_forward) {
    int min_conn = 0;
    return next_successor(seed, forward, to_nodes, to_forwards, to_forward, min_conn);
}

size_t Sorter::next_successor(node_idx_t seed, bool forward, vector<node_idx_t> &to_nodes, vector<bool> &to_forwards, bool &to_forward, int &min_conn) {
    auto remaining_nodes = [this](node_idx_t n, bool forward) {
        int count = 0;
        for(auto arc : (*nodes)[n].arcs)
            if(!(*visited)[arc.successor] && arc.forward == forward)
                count++;
        return count;
    };

    size_t best = 0;
    switch(extending_method){
        case extending_method_t::LESS_CONNECTED: {
            min_conn = INT32_MAX;
            for(size_t i = 0; i < to_nodes.size(); i++){
                int rn = remaining_nodes(to_nodes[i], to_forwards[i]);
                if(rn < min_conn) {
                    min_conn = rn;
                    best = i;
                }
            }
        }
            break;
        case extending_method_t::MORE_CONNECTED: {
            int max_conn = 0;
            for(size_t i = 0; i < to_nodes.size(); i++){
                int rn = remaining_nodes(to_nodes[i], to_forwards[i]);
                if(rn > max_conn) {
                    max_conn = rn;
                    best = i;
                }
            }
        }
            break;
        case extending_method_t::FIRST: // choose always the first
            // do nothing, it's before the cycle
            break;
        case extending_method_t::SIMILAR_ABUNDANCE: {
                uint32_t best_value = UINT32_MAX;
                for(size_t i = 0; i < to_nodes.size(); i++){
                    uint32_t ab_seed = nodes->at(seed).abundances.back();
                    uint32_t ab_succ = nodes->at(to_nodes.at(i)).abundances.front();

                    if(!forward)
                        ab_seed = nodes->at(seed).abundances.front();
                    if(!to_forwards.at(i))
                        ab_succ = nodes->at(to_nodes.at(i)).abundances.back();

                    // compute the distance
                    uint32_t diff = d(ab_seed, ab_succ);

                    if(diff == 0){ // same abundance!
                        best = i;
                        break;
                    }
                    if(diff < best_value){
                        best_value = diff;
                        best = i;
                    }
                }
            }
            break;
        case extending_method_t::SIMILAR_MEDIAN_ABUNDANCE:
            {
                auto best_value = UINT32_MAX;
                for(size_t i = 0; i < to_nodes.size(); i++){
                    auto ab_seed = nodes->at(seed).median_abundance;
                    auto ab_succ = nodes->at(to_nodes.at(i)).median_abundance;

                    // compute the distance
                    auto diff = d(ab_seed, ab_succ);

                    if(diff == 0){ // same abundance!
                        best = i;
                        break;
                    }
                    if(diff < best_value){
                        best_value = diff;
                        best = i;
                    }
                }
            }
            break;
        case extending_method_t::LOWER_MEDIAN_ABUNDANCE:{
                uint32_t min_ab = UINT32_MAX;
                for (size_t i = 0; i < to_nodes.size(); i++) {
                    auto ab = (*nodes)[to_nodes[i]].median_abundance;
                    if (ab < min_ab) {
                        min_ab = ab;
                        best = i;
                    }
                }
            }
            break;
        case extending_method_t::BIGGER_LENGTH:
            {
                uint32_t max_len = 0;
                for (size_t i = 0; i < to_nodes.size(); i++) {
                    auto len = (*nodes)[to_nodes[i]].length;
                    if (len > max_len) {
                        max_len = len;
                        best = i;
                    }
                }
            }
            break;
        case extending_method_t::SMALLER_LENGTH:
            {
                uint32_t min_len = UINT32_MAX;
                for (size_t i = 0; i < to_nodes.size(); i++) {
                    auto len = (*nodes)[to_nodes[i]].length;
                    if (len < min_len) {
                        min_len = len;
                        best = i;
                    }
                }
            }
            break;
        case extending_method_t::RANDOM:
            best = get_rand(to_nodes.size() - 1);
            break;
        default:
            cerr << "seed_successor(): unknown extending method!" << endl;
            exit(EXIT_FAILURE);
    }

    to_forward = to_forwards[best];
    last_node = to_nodes[best];
    return last_node;
}
