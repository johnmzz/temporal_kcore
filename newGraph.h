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
#include <numeric>
#include <sstream>

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
    struct Nbr {
        int v;
        int t;
        int eid;

        Nbr(int v, int t, int eid) : v(v), t(t), eid(eid) {}
    };

    // skyline
    struct Label {
        int v {};
        int start {};
        int end {};
        int act {};
        int last {};
        int prev {-1};
        int next {-1};
        bool aux {};

        Label(int s, int e) : start(s), end(e) {}

        // Label(const Label& other) {
        //     v = other.v;
        //     start = other.start;
        //     end = other.end;
        //     act = other.act;
        //     last = other.last;
        //     prev = other.prev;
        //     next = other.next;
        //     aux = other.aux;
        // }
    };

    unsigned int n_;    // num vertices
    unsigned int m_;    // num edges (temporal graph)
    int t_;             // number of time stamps
    int k_;
    int k_max_;
    int write_res_;
    unsigned long num_res_;
    unsigned long total_res_size_;
    unsigned long ct_size_;
    unsigned long skyline_size_;
    string version_;       // v or e

    int ts_;            // after conversion
    int te_;
    
    vector<pair<int, int>> edges_;          // edge list
    vector<vector<Nbr>> nbr_;    // graph adj list
    vector<int> edges_idx_;                 // index for timestamps of all edges

    vector<long> t_new_to_old_;
    unordered_map<long, int> t_old_to_new_;     // new t starts with 0

    vector<int> core_;      // cores of each v

    vector<unordered_map<int,int>> ct_cnt_;

    // core time
    vector<vector<pair<int,int>>> core_t_;
    vector<vector<pair<int,int>>> ct_e_;

    // skyline
    vector<vector<Label>> skyline_v_;
    vector<vector<Label>> skyline_;

    vector<int> last_;

    vector<vector<int>> ba_;
    vector<vector<int>> bs_;

    vector<Label> L;

    unordered_map<string, pair<int,int>> hres_;


public:
    void load(const string &path, long ts, long te, int k, string version, int write_res);
    void truncate(int ts, int te);

    // display
    void print_graph();
    void print_core();
    void print_ct();
    void print_ct_e();
    void print_skyline();
    void print_L();
    void print_ba();
    void print_bs();
    void print_ll();

    // inline
    bool invalid(int u);

    // compute core time
    void core_decomposition(int ts);
    void compute_core_deg(const int &t_s, vector<unordered_map<int, int>>& cd, vector<int>& core);
    void init_core_time(int ts, int te);
    void e_core_time();
    void core_time();
    void ct_size();
    void ct_e_size();

    //baseline
    void baseline();

    // distinct sets
    void ct_to_skyline();
    void init();
    void compute_sets();
    void enumerate(int h, int t1, int t2);
    void enumerate_e(int h, int t1, int t2);
    void output(int ts, int te, unordered_set<int>& R);

    // temporary (remove after use)
    void init_skyline();

    // linked list
    void remove(int x);
    void insert(int x, int a, int b);

};

inline bool Graph::invalid(int u) {
    if (core_[u] < k_) return true;
    return (!core_t_[u].empty()) && core_t_[u].back().second == te_+1;
}
#endif //TIME_RANGE_KCORE_H