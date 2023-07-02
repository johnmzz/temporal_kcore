#ifndef TIME_RANGE_KCORE_H
#define TIME_RANGE_KCORE_H

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
#include <iomanip>
#include <climits>

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
    int min_k_;
    int threads_;
    int k_max_{};
    int chunk_size_{};
    long t_min_{};
    long long idx_size_;
    unsigned int n_{};
    unsigned int m_{};
    unsigned int max_deg_{};
    unsigned int t_{};
    unsigned int effective_m_{};
    unsigned int max_effective_deg_{};

    bool* v_a_;
    bool* v_b_;
    int* t_offset_{};
    int* core_{};

    vector<int> ct_init_;
    vector<int> ctn_v2_;
    vector<long> t_new_to_old_;
    unordered_map<long, int> t_old_to_new_;
    vector<pair<int, int>> edges_;
    vector<int> edges_idx_;
    vector<int> offset_t_;
    vector<vector<pair<int, int>>> nbr_;
    unordered_map<int, int>* nbr_cnt_;
    unordered_map<int, int>* cd_;
    unordered_map<int, int>* ct_cnt_;
    vector<vector<pair<int,int>>>* core_t_{};
    // vector<int>* ct_;
    vector<pair<int,int>>* ct_;
    vector<vector<int>> old_ct_;
    vector<vector<int>> nbr_time_;
    vector<vector<int>> offset_;
    // vector<vector<unordered_map<int,int>>> ctn_;
    vector<unordered_set<int>>* ctn_;
    vector<vector<pair<int, int>>> messages_;


public:
    Graph();
    ~Graph();

    // construct graph
    void load(const string &graph_path);
    void init_nbr_cnt();
    void init_nbr_time();

    // inline
    bool invalid(int u, int k);
    bool invalid_v2(int u, int k, vector<pair<int,int>>* core_t);
    void del_nbr(int u, int v);
    void del_nbr_v2(int u, int v, vector<unordered_map<int,int>>& ct_cnt);

    // time range k-core
    void time_range_kcore(long _ts, long _te, int _k);
    void truncate(int ts, int te);
    void truncate_t(int ts, int te);
    void core_decomposition(int _k);
    void compute_core_deg(const int &t_s);
    void init_core_time(int ts, int te, int k);
    void init_ct_cnt(int k);

    // parallel ver
    // void time_range_kcore_parallel(long ts, long te, int _k, int threads);
    // void init_ct(int ts, int te, int k);
    // void init_cnt(int ts, int te, int k);
    // void local_ct(int u, int ts, int k, int &offset, vector<bool> &visited, vector<int> &nbr_t, vector<int> &bm_history);
    // // void local_ct(int u, int t_s, int t_e, int k);
    // void init_ctn(int ts, int te, int k);
    // int get_ct(int u, int ts);
    // void insert_ct(int u, int ts, int ct);

    // ver 2 (parallel init, incremental update)
    void time_range_kcore_v2(long _ts, long _te, int k, int threads);
    void init_core_time_v2(int _ts, int _te, int _k, int threads, vector<int>& offset, vector<unordered_map<int,int>>& ct_cnt, vector<pair<int,int>>* core_t);
    
    void init_ct_v2(int ts, int te, int k, vector<int> &offset, vector<int>& ct_init);
    void init_ctn_v2(int ts, int te, int k, vector<int> &offset, vector<int>& ct_init, vector<unordered_map<int,int>>& ct_cnt);
    void local_ct_v2(int u, int ts, int k, vector<bool> &visited, vector<int> &offset, vector<int>& ct_prev, vector<int>& ct_curr);

    // ver 3 (parallel naive, updated CTN of all neighbors each round)
    void time_range_kcore_v3(long _ts, long _te, int k, int threads);

    // display
    void print_ct(vector<pair<int,int>>* core_t);
    void print_ctn(vector<unordered_map<int,int>>& ct_cnt);
    void print_graph();
    void print_message();
    void print_nbr_cnt_();
    void print_idx_size();
    void print_local_ct();
    void print_nbr_time();
    void print_graph_size();
    void print_queue(queue<int> q);

    // output
    void write_index(vector<pair<int,int>>* core_t);
};

inline bool Graph::invalid(int u, int k) {
    if (core_[u] < k) return true;
    return (!core_t_[u][k].empty()) && core_t_[u][k].back().second == t_;
}

inline bool Graph::invalid_v2(int u, int k, vector<pair<int,int>>* core_t) {
    if (core_[u] < k) return true;
    return (!core_t[u].empty()) && core_t[u].back().second == t_;
}

inline void Graph::del_nbr(int u, int v) {
    if (ct_cnt_[u].find(v) == ct_cnt_[u].end()) return;
    --ct_cnt_[u][v];
    if (ct_cnt_[u][v]==0) ct_cnt_[u].erase(v);
}

inline void Graph::del_nbr_v2(int u, int v, vector<unordered_map<int,int>>& ct_cnt) {
    if (ct_cnt[u].find(v) == ct_cnt[u].end()) return;
    --ct_cnt[u][v];
    if (ct_cnt[u][v]==0) ct_cnt[u].erase(v);
}

#endif //TIME_RANGE_KCORE_H
