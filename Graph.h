//
// Created by DW on 7/29/2020.
//

#ifndef SPAN_CORE_GRAPH_H
#define SPAN_CORE_GRAPH_H


#include <string>
#include <vector>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <cstring>
#include <bitset>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <ctime>

#define _LINUX_

#ifdef _LINUX_
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#endif

using namespace std;

class Graph {

    unsigned int n_{};
    unsigned int m_{};
    unsigned int max_deg_{};
    unsigned int effective_m_{};
    int min_k_;
    long long idx_size_;

    unsigned int max_effective_deg_{};

    FILE* log_f_;

    vector<long> t_new_to_old_;
    vector<pair<int,int>> edges_;
    vector<int> edges_idx_;
    vector<vector<pair<int,int>>> nbr_;
//    vector<vector<int>> reverse_nbr_idx_;

    unordered_map<int,int>* nbr_cnt_;

    int* core_{};


    vector<int>* query_v_{};
    int* query_deg_{};

    unordered_map<int,int>* cd_;

    vector<vector<pair<int,int>>>* core_t_{};


    bool* v_a_;
    bool* v_b_;


    unordered_map<int,int>* ct_cnt_;

    int* t_offset_{};


    void init_nbr_cnt();
    void core_decomposition();
    void init_core_time();
    void init_ct_cnt(int k);

    void del_nbr(int u, int v);
    bool invalid(int u, int k);

//    functions for baseline index construction
    void decremental_core_bl(const int &t_s);
    void compute_core_time_bl(const int &t_s);
    void compute_core_deg(const int &t_s);

    void test_core_decomposition(const int &t_s);
    void print_idx_size();
    void print_graph_size();

public:
    unsigned int t_{};
    int k_max_{};
    Graph();
    ~Graph();
    void load(const string &path, bool timestamp_third);
    void load(const string &path, const int &total_edge_num, bool timestamp_third);

    bool query(int u, int t_s, int t_e, int k);
    int query_all(int t_s, int t_e, int k);
    void query_subgraph(int u, int t_s, int t_e, int k, vector<int>& r, vector<pair<int,int>>& r_edges);

    void query_init();

    void index();
    void index_baseline();


    void load_idx(const string &path);
    void write_idx(const string &path);
    void init_log(const string &log_path);

    void naive_index();
    void naive_index_size();

    void test();


//    void online_span_core(const int &u, const int &t_s, const int &t_e);
    void online_core_decomposition(const int &t_s, const int &t_e);
    int online_k_core(const int &t_s, const int &t_e, const int& k);
    int online_query(const int &t_s, const int &t_e, const int& k, int& snapshot_m);
    int online_query(const int &t_s, const int &t_e, const int& k);

//    for case study
    int online_span_core(const int &t_s, const int &t_e, const int& k);
    void edge_intersection(vector<vector<pair<int,int>>>& edges, vector<pair<int,int>>& result);
    int index_span_core(const int &t_s, const int &t_e, const int& k);
};




inline void Graph::del_nbr(int u, int v) {
//    if (core_[u] < k) return;
    if (ct_cnt_[u].find(v) == ct_cnt_[u].end()) return;
    --ct_cnt_[u][v];
    if (ct_cnt_[u][v]==0) ct_cnt_[u].erase(v);
}

inline bool Graph::invalid(int u, int k) {
    if (core_[u] < k) return true;
    return (!core_t_[u][k].empty()) && core_t_[u][k].back().second == t_;
}


bool cmp(const pair<int,int> &a, const pair<int,int> &b);
bool cmp_nbr(const pair<int,int> &a, const pair<int,int> &b);

#endif //SPAN_CORE_GRAPH_H
