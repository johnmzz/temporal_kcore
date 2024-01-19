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
#include <map>
#include <unordered_set>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <climits>
#include <chrono>

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
    unsigned int n_;    // num vertices
    unsigned int m_;    // num edges (temporal graph)
    int t_;             // number of time stamps
    int k_;
    int k_max_;
    // unsigned int effective_m_;
    // unsigned int max_effective_deg_;
    int num_cores_;     // number of distinct k-cores
    int threads_;
    int chunk_size_;    // divide t by chunks depend on number of threads
    int write_result_;
    int threads_init_;
    
    int ts_;            // after conversion
    int te_;

    // vector<int> core_;  // cores of each v
    vector<vector<int>> cores_;     // cores of each v for different threadss

    vector<vector<pair<int, int>>> nbr_;    // graph adj list
    vector<pair<int, int>> edges_;          // edge list
    vector<int> edges_idx_;                 // index for timestamps of all edges

    unordered_map<long, int> t_old_to_new_;     // new t starts with 0
    vector<long> t_new_to_old_;

    // vector<unordered_map<int, int>> nbr_cnt_;      // for each v: store {neighbor : count}

    vector<pair<pair<int,int>, vector<int>>> res_;  // store distinct k-cores

    vector<vector<pair<pair<int,int>, vector<int>>>> res_parallel_;     // store distinct k-cores parallel
    vector<int> num_cores_parallel_;

    vector<vector<pair<int, int>>> messages_;

    vector<vector<int>> offsets_;
    vector<vector<unordered_map<int,int>>> ct_cnts_;
    vector<vector<vector<pair<int,int>>>> core_ts_;

public:
    void load(const string &graph_path, long ts, long te, int k, int threads, int write_res);
    void truncate(int ts, int te);
    // void init_nbr_cnt(int ts);

    // display
    void print_graph();
    void print_ct(vector<vector<pair<int,int>>>& core_t);
    void print_ctn(vector<unordered_map<int, int>>& ct_cnt);

    // result (distinct k-cores)
    void t_cores(int ts, int te, int min_covered, vector<vector<pair<int,int>>>& core_t);
    void t_cores_parallel(int ts, int te, int min_covered, vector<vector<pair<int,int>>>& core_t, int part);

    // inline
    bool invalid(int u, vector<vector<pair<int,int>>>& core_t, int tid);

    // baseline serial (v1)
    void baseline_serial();
    void core_decomposition_org(int ts, int tid);
    void init_core_time(int ts, int te, vector<vector<pair<int,int>>>& core_t, int tid);
    void compute_core_deg(const int &t_s, vector<unordered_map<int, int>>& cd, vector<int>& core);

    // decrement parallel (v2)
    void decrement_parallel();

    // queue parallel (v3)
    void queue_parallel();

    // vertex-centric parallel (v4)
    void vcentric_parallel();
    void core_decomposition(int ts, int tid);
    void connect(int tid, int ts);
    void vcentric_init(int ts, vector<int>& offset, vector<unordered_map<int,int>>& ct_cnt, vector<vector<pair<int,int>>>& core_t, int partition_size, int tid);
    void init_ct(int ts, vector<int>& offset, vector<int>& ct_init, int tid);
    void init_ctn(int t_s, vector<int>& offset, vector<int>& ct_init, vector<unordered_map<int, int>>& ct_cnt, int tid);
    void local_ct(int u, int t_s, vector<bool>& visited, vector<int>& offset, vector<int>& ct_prev, vector<int>& ct_curr, int tid);

    // output
    void write_index(vector<vector<pair<int,int>>>& core_t);
    void write_result();
    void write_res_parallel();
};

inline bool Graph::invalid(int u, vector<vector<pair<int,int>>>& core_t, int tid) {
    if (cores_[tid][u] < k_) return true;
    return (!core_t[u].empty()) && core_t[u].back().second == te_+1;
}
#endif //TIME_RANGE_KCORE_H