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
    FILE* log_f_;

    int min_k_;
    int k_max_{};
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
    // vector<pair<int,int>>* ct_;
    vector<vector<bool>> v_c_;
    vector<vector<int>> nbr_time_;
    vector<vector<int>> offset_;


public:
    Graph();
    ~Graph();

    // construct graph
    void load(const string &graph_path);
    void init_nbr_cnt();
    void init_nbr_time();

    // inline
    bool invalid(int u, int k);
    void del_nbr(int u, int v);
    bool invalid_(int u, int k);
    int get_ct(int u, int ts);
    void insert(int u, int ts, int t);

    // time range k-core
    void time_range_kcore(long _ts, long _te, int _k);
    void truncate(int ts, int te);
    void core_decomposition(int _k);
    void compute_core_deg(const int &t_s);
    void init_core_time(int ts, int te, int k);
    void init_ct_cnt(int k);

    // parallel ver
    void time_range_kcore_parallel(long ts, long te, int _k);
    void local_ct(int u, int ts, int te, int k);

    // display
    void print_ct();
    void print_graph();
    void print_nbr_cnt_();
    void print_idx_size();
    void print_local_ct();
    void print_nbr_time();
    void print_graph_size();
    void print_queue(queue<int> q);
    void init_log(const string &log_path);
};

inline bool Graph::invalid(int u, int k) {
    if (core_[u] < k) return true;
    return (!core_t_[u][k].empty()) && core_t_[u][k].back().second == t_;
}

inline bool Graph::invalid_(int u, int k) {
    if (core_[u] < k) return true;
    return (!ct_[u].empty()) && ct_[u].back().second == t_;
}

inline void Graph::del_nbr(int u, int v) {
    if (ct_cnt_[u].find(v) == ct_cnt_[u].end()) return;
    --ct_cnt_[u][v];
    if (ct_cnt_[u][v]==0) ct_cnt_[u].erase(v);
}

inline int Graph::get_ct(int u, int ts) {
    for (int i = ct_[u].size()-1; i >= 0; --i) {
        if (ts >= ct_[u][i].first) {
            return ct_[u][i].second;
        }
    }
    return t_;
}

inline void Graph::insert(int u, int ts, int t) {
    auto idx = lower_bound(ct_[u].begin(), ct_[u].end(), make_pair(ts, 0));

    if (idx == ct_[u].end()) ct_[u].emplace_back(make_pair(ts, t));

    int i = idx - ct_[u].begin();
    if (ct_[u][i].first == u) {
        ct_[u][i].second = t;
    } else {
        ct_[u].insert(idx, make_pair(ts, t));
    }

    if (t == t_) {
        ct_[u].resize(i+1);
    }
}

#endif //TIME_RANGE_KCORE_H
