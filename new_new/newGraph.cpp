#include "newGraph.h"
#include <omp.h>

void Graph::print_graph() {
    cout << "nbr_ = " << endl;
    for (int u = 0; u < nbr_.size(); u++) {
        if (!nbr_[u].empty()) {
            cout << u << ": [";
            for (auto vt : nbr_[u]) {
                cout << "<" << vt.first << "," << vt.second << ">,";
            }            
            cout << "]" << endl;
        }
    }
}

void Graph::print_ct(vector<vector<pair<int,int>>>& core_t) {
    cout << "core_t = " << endl;
    for (int u = 0; u < n_; u++) {
        cout << "u = " << u << ": [";
        for (auto v : core_t[u]) {
            cout << "<" << v.first << "," << v.second << ">,";
        }
        cout << "]," << endl;
    }
}

void Graph::print_ctn(vector<unordered_map<int, int>>& ct_cnt) {
    cout << "ct_cnt = " << endl;
    for (int u = 1; u < n_; u++) {
        cout << std::setw(2) << "u = " << u << ": [";
        if (ct_cnt[u].empty()) {
            cout << "]\n";
            continue;
        }
        for (auto it = ct_cnt[u].begin(); it != ct_cnt[u].end(); it++) {
            cout << "<" << it->first << "," << it->second << ">,";
        }
        cout << "\n";
    }
}

void Graph::write_index(vector<vector<pair<int,int>>>& core_t) {
    cout << "writing to 1.txt\n";
    ofstream fout("1.txt");

    for (int u = 0; u < n_; u++) {
        fout << "u = " << u << ": [";
        for (auto v : core_t[u]) {
            fout << "<" << v.first << "," << v.second << ">,";
        }
        fout << "]" << endl;
    }

    fout.close();
}

void Graph::write_result() {
    cout << "writing result to res.txt\n";
    ofstream fout("res.txt");

    for (auto p : res_) {
        auto time_pair = p.first;
        auto tcore = p.second;
        fout << "<" << time_pair.first << "," << time_pair.second << ">: [";
        
        for (auto v : tcore) {
            fout << v << ",";
        }
        fout << "]\n";
    }

    fout.close();
}

void Graph::load(const string &path, long _ts, long _te, int k, int threads, int write_res) {
    printf("Graph path: %s\n", path.c_str());

    if (write_res) {
        write_result_ = true;
    } else {
        write_result_ = false;
    }

    ifstream ifs(path);
    if (!ifs.is_open()) {
        cerr << "open file failed!" << endl;
        exit(1);
    }

    n_ = 0;
    m_ = 0;
    k_ = k;
    threads_ = threads;
    threads_init_ = 1;

    if (path == "../data/youtube_sorted.txt" || path == "../data/wikipedia.txt") {
        threads_init_ = threads;
        if (threads > 4) threads_init_ = 4;
    }
    
    int u,v;
    long ts;
    while (ifs.good() && !ifs.eof()) {
        char line[200];
        ifs.getline(line, 200);
        if (line[0] < '0' || line[0] > '9') {
            continue;
        }
        sscanf(line, "%d %d %ld", &u, &v, &ts);

        if (u == v) {
            continue;
        }

        edges_.emplace_back(make_pair(u,v));

        // adjust size of neighbor list if necessary
        if (u+1 > nbr_.size()) {
            nbr_.resize(u+1);
        }
        if (v+1 > nbr_.size()) {
            nbr_.resize(v+1);
        }

        long pre_ts = -1;
        if (!t_new_to_old_.empty()) {
            pre_ts = t_new_to_old_.back();
        }

        if (ts != pre_ts) {
            t_new_to_old_.emplace_back(ts);
            t_old_to_new_.insert(make_pair(ts, t_new_to_old_.size()-1));
            edges_idx_.emplace_back(edges_.size() - 1);
        }

        int format_t = t_new_to_old_.size() - 1;

        // construct neighbor list
        nbr_[u].emplace_back(make_pair(v,format_t));
        nbr_[v].emplace_back(make_pair(u,format_t));
    }
    ifs.close();

    n_ = nbr_.size();
    edges_idx_.emplace_back(edges_.size());

    ts_ = t_old_to_new_[_ts];
    te_ = t_old_to_new_[_te];
    truncate(ts_, te_);

    // init_nbr_cnt(ts_);
}

void Graph::truncate(int ts, int te) {
    t_ = te - ts + 1;
    vector<vector<pair<int, int>>> nbr(n_);
    for (long t = ts; t <= te; t++) {
        for (int i = edges_idx_[t]; i < edges_idx_[t+1]; i++) {
            int u = edges_[i].first;
            int v = edges_[i].second;

            nbr[u].emplace_back(make_pair(v,t));
            nbr[v].emplace_back(make_pair(u,t));
        }
    }
    nbr_ = nbr;
}

void Graph::t_cores(int ts, int te, int min_covered, vector<vector<pair<int,int>>>& core_t) {
    map<int, vector<int>> H;    // TODO: optimize
    for (int u = 0; u < n_; u++) {
        if (!core_t[u].empty() && core_t[u].back().second != te+1) {
            H[core_t[u].back().second].emplace_back(u);
        }
    }

    vector<int> t_core;
    for (auto it = H.lower_bound(min_covered); it != H.end(); it++) {
        t_core.insert(t_core.end(), it->second.begin(), it->second.end());
        num_cores_++;
        if (write_result_) res_.emplace_back(make_pair(make_pair(ts, it->first), t_core));
    }
}

void Graph::t_cores_parallel(int ts, int te, int min_covered, vector<vector<pair<int,int>>>& core_t, int part) {
    map<int, vector<int>> H;
    for (int u = 0; u < n_; u++) {
        if (!core_t[u].empty() && core_t[u].back().second != te+1) {
            H[core_t[u].back().second].emplace_back(u);
        }
    }

    vector<int> t_core;
    for (auto it = H.lower_bound(min_covered); it != H.end(); it++) {

        t_core.insert(t_core.end(), it->second.begin(), it->second.end());
        num_cores_parallel_[part]++;
        // if (write_result_) res_.emplace_back(make_pair(make_pair(ts, it->first), t_core));
        if (write_result_) res_parallel_[part].emplace_back(make_pair( make_pair(ts, it->first), t_core));;
    }
}

void Graph::core_decomposition_org(int ts, int tid) {
    // init_nbr_cnt
    int effective_m = 0;
    int max_effective_deg = 0;

    vector<unordered_map<int,int>> nbr_cnt(n_);

    for (int u = 0; u < n_; ++u) {
        for(auto &i : nbr_[u]) {
            if (i.second < ts) continue;
            if (nbr_cnt[u].find(i.first) != nbr_cnt[u].end()) {
                ++nbr_cnt[u][i.first];
            } else {
                nbr_cnt[u].insert(make_pair(i.first,1));
            }
        }
        if(nbr_cnt[u].size() > max_effective_deg) max_effective_deg = nbr_cnt[u].size();
        effective_m += nbr_cnt[u].size();
    }
    effective_m /= 2;

    vector<int> core(n_);
    
    vector<int> t_offset(n_);
    vector<int> vert(n_);
    vector<int> bin(max_effective_deg+1);

    vector<bool> v_a(n_, false);

    for (int u = 0; u < n_; ++u) {
        int d = nbr_cnt[u].size();
        core[u] = d;
        ++bin[d];
    }

    // sort vertices based on degree
    int offset = 0;
    for (int i = 0; i <= max_effective_deg ; ++i) {
        int num = bin[i];
        bin[i] = offset;
        offset += num;
    }
    for (int u = 0; u < n_; ++u) {
        t_offset[u] = bin[core[u]];
        vert[t_offset[u]] = u;
        bin[core[u]]++;
    }

    for (int i = max_effective_deg; i >= 1; --i) bin[i] = bin[i - 1];
    bin[0] = 0;
    k_max_ = 0;

    for (int i = 0; i < n_; ++i) {
        int u = vert[i];

        for (auto& item : nbr_[u]){
            if (v_a[item.first] || item.second < ts) continue;
            v_a[item.first] = true;
            if (core[item.first] > core[u]){
                int dv = core[item.first], pv = t_offset[item.first];
                int pw = bin[dv], w = vert[pw];
                if (item.first != w){
                    t_offset[item.first] = pw, vert[pv] = w;
                    t_offset[w] = pv, vert[pw] = item.first;
                }
                ++bin[dv];
                --core[item.first];
            }
        }

        for (auto& item : nbr_[u]){
            v_a[item.first] = false;
        }

        if (core[u] > k_max_) k_max_ = core[u];     // TODO: atomic?
    }

    cores_[tid] = core;     // TODO: atomic?
}

void Graph::core_decomposition(int ts, int tid) {
    // init_nbr_cnt
    int effective_m = 0;
    int max_effective_deg = 0;

    vector<unordered_map<int,int>> nbr_cnt(n_);

    for (int u = 0; u < n_; ++u) {
        for(auto &i : nbr_[u]) {
            if (i.second < ts) continue;
            if (nbr_cnt[u].find(i.first) != nbr_cnt[u].end()) {
                ++nbr_cnt[u][i.first];
            } else {
                nbr_cnt[u].insert(make_pair(i.first,1));
            }
        }
        if(nbr_cnt[u].size() > max_effective_deg) max_effective_deg = nbr_cnt[u].size();
        effective_m += nbr_cnt[u].size();
    }
    effective_m /= 2;

    vector<int> core(n_);
    vector<int> t_offset(n_);
    vector<int> vert(n_);
    vector<int> bin(max_effective_deg+1);

    vector<bool> v_a(n_, false);

    for (int u = 0; u < n_; ++u) {
        int d = nbr_cnt[u].size();
        core[u] = d;
        ++bin[d];
    }

    // sort vertices based on degree
    int offset = 0;
    for (int i = 0; i <= max_effective_deg ; ++i) {
        int num = bin[i];
        bin[i] = offset;
        offset += num;
    }
    for (int u = 0; u < n_; ++u) {
        t_offset[u] = bin[core[u]];
        vert[t_offset[u]] = u;
        bin[core[u]]++;
    }

    for (int i = max_effective_deg; i >= 1; --i) bin[i] = bin[i - 1];
    bin[0] = 0;
    k_max_ = 0;

    int i;
    for (i = 0; i < n_; ++i) {
        int u = vert[i];
        if (core[u] >= k_) {
            k_max_ = core[u];
            break;
        }

        for (auto& item : nbr_[u]){
            if (v_a[item.first] || item.second < ts) continue;
            v_a[item.first] = true;
            if (core[item.first] > core[u]){
                int dv = core[item.first], pv = t_offset[item.first];
                int pw = bin[dv], w = vert[pw];
                if (item.first != w){
                    t_offset[item.first] = pw, vert[pv] = w;
                    t_offset[w] = pv, vert[pw] = item.first;
                }
                ++bin[dv];
                --core[item.first];
            }
        }

        for (auto& item : nbr_[u]){
            v_a[item.first] = false;
        }

        if (core[u] > k_max_) k_max_ = core[u];     // TODO: atomic?
    }
    for (;i < n_; i++) {
        int u = vert[i];
        core[u] = k_;
    }

    cores_[tid] = core;     // TODO: atomic?
}

void Graph::compute_core_deg(const int &t_s, vector<unordered_map<int, int>>& cd, vector<int>& core) {
    for (int u = 0; u < n_; ++u) {
        cd[u].clear();
        for (int i = nbr_[u].size()-1;i>=0;--i){
            int v = nbr_[u][i].first;
            int t = nbr_[u][i].second;
            if (t < t_s) continue;

            if (core[v] < core[u]) continue;

            if (cd[u].find(v) == cd[u].end()){
                cd[u].insert(make_pair(v,1));
            }else{
                ++cd[u][v];
            }
        }
    }
}

void Graph::init_core_time(int _ts, int _te, vector<vector<pair<int,int>>>& core_t, int tid) {
    vector<int> core = cores_[tid];

    vector<unordered_map<int,int>> cd(n_);
    compute_core_deg(_ts, cd, core);

    vector<bool> v_a(n_, false);
    vector<bool> v_b(n_, false);
    queue<int> q;
    vector<int> cnt(k_max_+1);
    for (int t_e = _te; t_e >= _ts; --t_e) {
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e+1]; ++i) { // delete edge
            int u = edges_[i].first;
            int v = edges_[i].second;

            // process u
            if (core[u] <= core[v]){
                --cd[u][v];
                if (cd[u][v] == 0){
                    cd[u].erase(v);
                    if (cd[u].size() < core[u] && !v_a[u]){
                        q.push(u);
                        v_a[u] = true;
                    }
                }
            }

            // process v
            if (core[v] <= core[u]){
                --cd[v][u];
                if (cd[v][u] == 0){
                    cd[v].erase(u);
                    if (cd[v].size() < core[v] && !v_a[v]){
                        q.push(v);
                        v_a[v] = true;
                    }
                }
            }
        }
        while (!q.empty()){
            int u = q.front();
            q.pop();
            v_a[u] = false;

            int oc = core[u];
            fill(cnt.begin(), cnt.end(), 0);

            // LocalCore
            for (int i = 0; i < nbr_[u].size(); ++i) {
                int v = nbr_[u][i].first;
                int t = nbr_[u][i].second;
                if (t >= t_e) break;
                if (t < _ts) continue;

                if (v_b[v]) continue;
                v_b[v] = true;

                ++cnt[core[v] < core[u] ? core[v]:core[u]];
            }
            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b[nbr_[u][i].first] = false;
            }

            int _cd = 0;
            for (int k = oc; k >= 0; --k) {
                _cd += cnt[k];
                if(_cd >= k){
                    core[u] = k;
                    break;
                }
            }

            // update cd;
            cd[u].clear();
            for (int i = 0; i < nbr_[u].size(); ++i) {
                int v = nbr_[u][i].first;
                int t = nbr_[u][i].second;
                if (t >= t_e) break;
                if (t < _ts) continue;
                if (core[v] < core[u]) continue;
                if (cd[u].find(v) == cd[u].end()) cd[u].insert(make_pair(v,1));
                else ++cd[u][v];
            }

            // add affected neighbor to the queue
            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                if (nbr_[u][i].second < _ts) continue;
                int v = nbr_[u][i].first;
                if (core[u] < core[v] && core[v] <= oc && !v_b[v]){    // v affected if oldcore(u) > core(v) but now core(u) < core(v)
                    v_b[v] = true;
                    cd[v].erase(u);
                    if (!v_a[v] && cd[v].size() < core[v]){
                        q.push(v);
                        v_a[v] = true;
                    }
                }
            }
            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b[nbr_[u][i].first] = false;
            }

            for (int k = oc; k > core[u]; --k) {
                if (k == k_) {
                    core_t[u].emplace_back(make_pair(_ts, t_e)); // if u popped from q, then its core changed at t_e, therefore put in Index
                }
            }
        }
    }
}

void Graph::init_ct(int ts, vector<int>& offset, vector<int>& ct_init, int tid) {
    vector<int> bm_history;
    vector<bool> visited(n_, false);

    for (int u = 1; u < n_; ++u) {
        if (cores_[tid][u] < k_) continue;

        int cnt = 0;
        for (int i = 0; i < nbr_[u].size(); ++i) {
            int t = nbr_[u][i].second;
            if (t < ts) {
                offset[u] = i+1;    // thread local variable
                continue;
            }

            int v = nbr_[u][i].first;
            if (visited[v]) continue;   // function local variable

            visited[v] = true;
            bm_history.emplace_back(v);

            ++cnt;
            if (cnt == k_) {
                ct_init[u] = t;
                break;
            }
        }
        for (auto &v : bm_history) visited[v] = false;
        bm_history.clear();
    }
}

void Graph::init_ctn(int t_s, vector<int>& offset, vector<int>& ct_init, vector<unordered_map<int, int>>& ct_cnt, int tid) {
    for (int u = 1; u < n_; ++u) {
        if (cores_[tid][u] < k_) continue;
        
        int t = ct_init[u];
        for (int i = offset[u]; i < nbr_[u].size(); ++i){
            // if (nbr_[u][i].second < t_s) continue;
            if (nbr_[u][i].second > t) break;

            int v = nbr_[u][i].first;
            if (cores_[tid][v] < k_ || t < ct_init[v]) continue;   // v's CT < u's CT, then v is in k-core earlier than u, thus add
            
            if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                ct_cnt[u].insert(make_pair(v,1));
            }else{
                ++ct_cnt[u][v];
            }
        }
    }
}

void Graph::local_ct(int u, int t_s, vector<bool>& visited, vector<int>& offset, vector<int>& ct_prev, vector<int>& ct_curr, int tid) {
    vector<int> nbr_t;
    vector<int> bm_history;
    int ub = 0;

    for (int i = offset[u]; i < nbr_[u].size(); ++i) {
        int t = nbr_[u][i].second;

        if (nbr_t.size() >= k_ && t > ub) break;

        // if (t < t_s) continue;

        int v = nbr_[u][i].first;
        if (cores_[tid][v] < k_ || visited[v] || ct_prev[v] == te_+1) continue;

        visited[v] = true;
        int v_t = ct_prev[v];
        int ct = max(t, v_t);

        nbr_t.emplace_back(ct);
        bm_history.emplace_back(v);

        if (nbr_t.size() <= k_) ub = max(ub, ct);
    }
    if (nbr_t.size() >= k_) {
        nth_element(nbr_t.begin(),nbr_t.begin()+k_-1,nbr_t.end());
        ct_curr[u] = nbr_t[k_-1];
    } else {
        ct_curr[u] = te_+1;
    }

    for (auto &v : bm_history) visited[v] = false;
}

void Graph::connect(int tid, int ts) {
    int partition_size;
    if (n_ % threads_ == 0) {
        partition_size = n_ / threads_;
    } else {
        partition_size = (n_ + (threads_ - n_ % threads_)) / threads_;
    }

    // declare variables
    vector<int> offset(n_);
    vector<unordered_map<int,int>> ct_cnt(n_);
    vector<vector<pair<int,int>>> core_t(n_);
    messages_.resize(threads_);

    core_decomposition(ts, tid);

    vcentric_init(ts, offset, ct_cnt, core_t, partition_size, tid);

    offsets_[tid] = offset;
    ct_cnts_[tid] = ct_cnt;
    core_ts_[tid] = core_t;
}

void Graph::vcentric_init(int ts, vector<int>& offset, vector<unordered_map<int,int>>& ct_cnt, vector<vector<pair<int,int>>>& core_t, int partition_size, int tid) {
    vector<int> ct_old(n_, te_+1);

    init_ct(ts, offset, ct_old, tid);
    init_ctn(ts, offset, ct_old, ct_cnt, tid);

    vector<int> ct_new = ct_old;

    int round = 0;
    bool update = true;
    while (update) {
        update = false;
        round++;

        #pragma omp parallel num_threads(threads_)
        {
            #pragma omp for schedule(static, partition_size)
            for (int u = 0; u < n_; ++u) {
                if (cores_[tid][u] < k_ || ct_old[u] == te_+1 || ct_cnt[u].size() >= k_) continue;

                vector<bool> visited(n_, false);

                int old_ct = ct_old[u];
                local_ct(u, ts, visited, offset, ct_old, ct_new, tid);
                int new_ct = ct_new[u];

                // re-compute ct_cnt[u]
                ct_cnt[u].clear();
                for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                    int t = nbr_[u][i].second;
                    // if (t < ts) continue;
                    if (t > new_ct) break;

                    int v = nbr_[u][i].first;
                    if (cores_[tid][v] < k_) continue;

                    if (ct_old[v] > new_ct) continue;

                    if (new_ct != te_+1) {
                        if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                            ct_cnt[u].insert(make_pair(v,1));
                        }else{
                            ++ct_cnt[u][v];
                        }
                    }

                    // update neighbor
                    if (cores_[tid][v] < k_ || visited[v]) continue;
                    if (new_ct > ct_old[v]) {
                        visited[v] = true;
                        #pragma omp critical
                        {
                            messages_[v/partition_size].emplace_back(make_pair(v,u));
                        }
                    }
                }
            }

            for (auto e : messages_[omp_get_thread_num()]) {
                int u = e.first;
                int v = e.second;
                //   received > self
                if (ct_new[v] > ct_new[u]) {    // if new CT[v] is later than CT[u] now, then v is not CT neighbor of u anymore
                    if (ct_cnt[u].find(v) != ct_cnt[u].end()) {
                        ct_cnt[u].erase(v);
                    }
                }
            }
            messages_[omp_get_thread_num()].clear();

            #pragma omp barrier

            #pragma omp for schedule(static, partition_size)
            for (int u = 0; u < n_; ++u) {
                if (ct_old[u] != ct_new[u]) {
                    ct_old[u] = ct_new[u];
                    update = true;
                }
            }
        }
    }
    for (int u = 1; u < n_; ++u) {
        if (cores_[tid][u] >= k_) {
            core_t[u].emplace_back(make_pair(ts, ct_new[u]));
        }
    }
    printf("Initialization complete, round taken = %d\n", round);
}







// v1
void Graph::baseline_serial() {
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    int tid = 0;
    cores_.resize(1);

    core_decomposition_org(ts_, tid);
    if (k_max_ < k_) {
        printf("queried k = %d exceed maximum core in G[ts,te] = %d\n", k_, k_max_);
        return;
    }

    // declare variables
    vector<bool> v_a(n_, false);
    vector<bool> v_b(n_, false);
    vector<int> offset(n_);

    vector<unordered_map<int,int>> ct_cnt(n_);
    vector<vector<pair<int,int>>> core_t(n_);

    // start-anchored core time
    init_core_time(ts_, te_, core_t, tid);

    // compute distincr k-cores
    t_cores(ts_, te_, -1, core_t);

    // compute ct_cnt before start
    for (int u = 0; u < n_; u++) {
        if (invalid(u,core_t,tid)) continue;
        ct_cnt[u].clear();

        int t = core_t[u].front().second;   // the new core time of u
        for (int i = offset[u]; i < nbr_[u].size(); ++i){
            if (nbr_[u][i].second> t) break;
            if (nbr_[u][i].second < ts_) continue;

            int v = nbr_[u][i].first;
            if (cores_[tid][v]<k_ || t < core_t[v].front().second) continue;
            if (ct_cnt[u].find(v) == ct_cnt[u].end()) {
                ct_cnt[u].insert(make_pair(v,1));
            }else{
                ++ct_cnt[u][v];
            }
        }
    }

    queue<int> q;
    bool update_flag = false;
    for (int t_s = ts_+1; t_s <= te_; ++t_s) {
        int min_covered = te_+1;
        for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
            int u = edges_[i].first;
            int v = edges_[i].second;

            if (invalid(u,core_t,tid) || invalid(v,core_t,tid)) continue;
            if (!v_a[u]){   // process u
                if (ct_cnt[u].find(v) != ct_cnt[u].end()) {
                    --ct_cnt[u][v];
                    if (ct_cnt[u][v] == 0) ct_cnt[u].erase(v);
                }
                if (ct_cnt[u].size() < k_){
                    q.push(u);
                    v_a[u] = true;
                }
            }
            if (!v_a[v]) {  // process v
                if (ct_cnt[v].find(u) != ct_cnt[v].end()) {
                    --ct_cnt[v][u];
                    if (ct_cnt[v][u] == 0) ct_cnt[v].erase(u);
                }
                if (ct_cnt[v].size() < k_) {
                    q.push(v);
                    v_a[v] = true;
                }
            }
        }
        if (!q.empty()) update_flag = true;

        while (!q.empty()){
            int u = q.front();
            q.pop();
            v_a[u] = false;

            // LocalCT
            vector<int> nbr_t;
            vector<int> bm_history;
            int ct = 0;
            for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                int t = nbr_[u][i].second;
                if (nbr_t.size() >= k_ && t > ct) break;
                if (t < t_s){   // ts_ must be less than t_s
                    offset[u] = i+1;
                    continue;
                }
               int v = nbr_[u][i].first;
                if (invalid(v,core_t,tid) || v_b[v]) continue;
                v_b[v] = true;
                int v_t = core_t[v].back().second;
                nbr_t.emplace_back(max(t,v_t));
                bm_history.emplace_back(v);

                if (nbr_t.size() <= k_) ct = max(ct,v_t);
            }
            for (auto &v:bm_history) v_b[v] = false;
            int new_t = te_+1;
            if (nbr_t.size() >= k_){
                nth_element(nbr_t.begin(),nbr_t.begin()+k_-1,nbr_t.end());
                new_t = nbr_t[k_-1];
            }

            // insert into index
            int old_t = core_t[u].back().second;
            min_covered = min(old_t, min_covered);
            if (core_t[u].back().first == t_s){
                core_t[u].back().second = new_t;
            }else{
                core_t[u].emplace_back(make_pair(t_s,new_t));
            }

            // re-compute ct_cnt[u]
            ct_cnt[u].clear();
            for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                int t = nbr_[u][i].second;
                int v = nbr_[u][i].first;
                if (t > new_t) break;
                if (t < ts_) continue;  // should never be triggered as i starts at offset[u], which should be set to position > ts_
                if (invalid(v,core_t,tid) || core_t[v].back().second > new_t) continue;

                if (new_t != te_+1){     // if ct(u) invalid, no need to add into index, but need to update neighbors
                    if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                        ct_cnt[u].insert(make_pair(v,1));
                    }else{
                        ++ct_cnt[u][v];
                    }
                }

                // add neighbor to queue if necessary
                if (v_a[v]) continue;
                if (core_t[v].back().second < old_t || new_t <= core_t[v].back().second) continue;
                ct_cnt[v].erase(u);
                if (ct_cnt[v].size() < k_){
                    q.push(v);
                    v_a[v] = true;
                }
            }
        }

        // compute distinct k-cores
        if (update_flag) {
            t_cores(t_s, te_, min_covered, core_t);
            update_flag = false;
        }
    }
    #ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
        printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);

        struct rusage rUsage;
        getrusage(RUSAGE_SELF, &rUsage);
        long ms = rUsage.ru_maxrss;
        printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
    #endif

    // write_index(core_t);
    cout << "number of results = " << num_cores_ << endl;
    // write_result();
}






// v2
void Graph::decrement_parallel(){
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    // calculate chunk size
    if (t_ % threads_ == 0) {
        chunk_size_ = t_ / threads_;
    } else {
        chunk_size_ = (t_ + (threads_ - t_ % threads_)) / threads_;
    }

    // set start times
    vector<int> start_time;
    for (int i = 0; i * chunk_size_ < t_; ++i) {
        start_time.emplace_back(ts_ + i * chunk_size_);
    }

    // variables
    cores_.resize(threads_);
    core_ts_.resize(start_time.size());
    res_parallel_.resize(start_time.size());
    num_cores_parallel_.resize(start_time.size());

    // init core times with different threads
    #pragma omp parallel for num_threads(threads_init_) schedule(static)
    for (int tid = 0; tid < start_time.size(); tid++) {
        vector<vector<pair<int,int>>> core_t(n_);

        // start-anchored core time
        core_decomposition_org(start_time[tid], tid);
        init_core_time(start_time[tid], te_, core_t, tid);
        
        core_ts_[tid] = core_t;  // TODO: atomic?

        t_cores_parallel(start_time[tid], te_, -1, core_ts_[tid], tid);
    }

    // final result
    vector<vector<pair<int,int>>> core_t_final(n_);

    // MAIN PROCESS
    #pragma omp parallel for num_threads(threads_) schedule(static) ordered
    for (int stidx = 0; stidx < start_time.size(); stidx++) {
        // allocate start and end time for each thread
        int tid = omp_get_thread_num();
        int t_start = start_time[stidx];
        int t_end;
        if (stidx == start_time.size()-1) {
            t_end = te_+1;
        } else {
            t_end = start_time[stidx+1];
        }

        // declare variables
        vector<bool> v_a(n_, false);
        vector<bool> v_b(n_, false);

        vector<int> offset(n_);
        vector<unordered_map<int,int>> ct_cnt(n_);
        vector<vector<pair<int,int>>> core_t = core_ts_[stidx];

        // re-compute ct_cnt before start
        for (int u = 0; u < n_; ++u) {
            if (invalid(u, core_t, tid)) continue;
            ct_cnt[u].clear();
            
            int t = core_t[u].front().second;
            for (int i = offset[u]; i < nbr_[u].size(); ++i){
                if (nbr_[u][i].second > t) break;
                if (nbr_[u][i].second < t_start) continue;

                int v = nbr_[u][i].first;
                if (cores_[tid][v] < k_ || core_t[v].empty() || t < core_t[v].front().second) continue;
                if (ct_cnt[u].find(v) == ct_cnt[u].end()) {
                    ct_cnt[u].insert(make_pair(v,1));
                }else{
                    ++ct_cnt[u][v];
                }
            }
        }
        
        queue<int> q;
        bool update_flag = false;
        for (int t_s = t_start+1; t_s < t_end; ++t_s) {
            int min_covered = te_+1;
            for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
                int u = edges_[i].first;
                int v = edges_[i].second;

                if (invalid(u,core_t,tid) || invalid(v,core_t,tid)) continue;
                if (!v_a[u]){   // process u
                    if (ct_cnt[u].find(v) != ct_cnt[u].end()) {
                        --ct_cnt[u][v];
                        if (ct_cnt[u][v] == 0) ct_cnt[u].erase(v);
                    }
                    if (ct_cnt[u].size()<k_){
                        q.push(u);
                        v_a[u] = true;
                    }
                }
                if (!v_a[v]) {  // process v
                    if (ct_cnt[v].find(u) != ct_cnt[v].end()) {
                        --ct_cnt[v][u];
                        if (ct_cnt[v][u] == 0) ct_cnt[v].erase(u);
                    }
                    if (ct_cnt[v].size() < k_) {
                        q.push(v);
                        v_a[v] = true;
                    }
                }
            }
            if (!q.empty()) update_flag = true;

            // process queue
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                v_a[u] = false;

                // LocalCT
                vector<int> nbr_t;
                vector<int> bm_history;
                int ct = 0;
                for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                    int t = nbr_[u][i].second;
                    if (nbr_t.size() >= k_ && t > ct) break;
                    if (t < t_s){
                        offset[u] = i+1;
                        continue;
                    }

                    int v = nbr_[u][i].first;
                    if (invalid(v,core_t,tid) || v_b[v]) continue;
                    v_b[v] = true;
                    int v_t = core_t[v].back().second;
                    nbr_t.emplace_back(max(t,v_t));
                    bm_history.emplace_back(v);

                    if (nbr_t.size() <= k_) ct = max(ct,v_t);
                }
                for (auto &v:bm_history) v_b[v] = false;
                int new_t = te_+1;
                if (nbr_t.size() >= k_){
                    nth_element(nbr_t.begin(),nbr_t.begin()+k_-1,nbr_t.end());
                    new_t = nbr_t[k_-1];
                }

                // insert into index
                int old_t = core_t[u].back().second;
                min_covered = min(old_t, min_covered);
                if (core_t[u].back().first == t_s){
                    core_t[u].back().second = new_t;
                }else{
                    core_t[u].emplace_back(make_pair(t_s,new_t));
                }

                // re-compute ct_cnt[u]
                ct_cnt[u].clear();
                for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                    int t = nbr_[u][i].second;
                    int v = nbr_[u][i].first;
                    if (t > new_t) break;
                    if (t < t_start) continue; // should never reach
                    if (invalid(v,core_t,tid) || core_t[v].back().second > new_t) continue;

                    if (new_t != te_+1){
                        if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                            ct_cnt[u].insert(make_pair(v,1));
                        }else{
                            ++ct_cnt[u][v];
                        }
                    }

                    // add neighbor to queue if necessary
                    if (v_a[v]) continue;
                    if (core_t[v].back().second < old_t || new_t <= core_t[v].back().second) continue;
                    ct_cnt[v].erase(u);
                    if (ct_cnt[v].size() < k_){
                        q.push(v);
                        v_a[v] = true;
                    }
                }
            }
            // compute distinct k-cores
            if (update_flag) {
                t_cores_parallel(t_s, te_, min_covered, core_t, stidx);
                update_flag = false;
            }
        }
        #pragma omp ordered
        {
            for (int u = 0; u < n_; u++) {
                for (int i = 0; i < core_t[u].size(); i++) {
                    if (!core_t_final[u].empty() && i == 0 && core_t_final[u].back().second == core_t[u][i].second) {
                        continue;
                    } else {
                        core_t_final[u].emplace_back(core_t[u][i]);
                    }
                }
            }
        }
    }

    #ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
        printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);

        struct rusage rUsage;
        getrusage(RUSAGE_SELF, &rUsage);
        long ms = rUsage.ru_maxrss;
        printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
    #endif

    // write_index(core_t_final);
    // write_res_parallel();
    // print_res_size();
    int total_num_cores = 0;
    for (auto v : num_cores_parallel_) total_num_cores += v;
    cout << "number of results = " << total_num_cores << endl;
}






// v3
void Graph::queue_parallel() {
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    // calculate chunk sizes
    if (n_ % threads_ == 0) {
        chunk_size_ = n_ / threads_;
    } else {
        chunk_size_ = (n_ + (threads_ - n_ % threads_)) / threads_;
    }

    int init_tid = 0;
    cores_.resize(1);

    core_decomposition_org(ts_, init_tid);
    if (k_max_ < k_) {
        printf("queried k = %d exceed maximum core in G[ts,te] = %d\n", k_, k_max_);
        return;
    }

    // declare variables
    vector<bool> v_a(n_, false);
    vector<int> offset(n_);
    vector<unordered_map<int,int>> ct_cnt(n_);
    vector<vector<pair<int,int>>> core_t(n_);
    messages_.resize(threads_);

    // start-anchored core time
    init_core_time(ts_, te_, core_t, init_tid);

    // compute distincr k-cores
    t_cores(ts_, te_, -1, core_t);

    // compute ct_cnt before start
    for (int u = 0; u < n_; u++) {
        if (invalid(u,core_t,init_tid)) continue;
        ct_cnt[u].clear();

        int t = core_t[u].front().second;   // the new core time of u
        for (int i = offset[u]; i < nbr_[u].size(); ++i){
            if (nbr_[u][i].second> t) break;
            if (nbr_[u][i].second < ts_) continue;

            int v = nbr_[u][i].first;
            if (cores_[init_tid][v]<k_ || t < core_t[v].front().second) continue;
            if (ct_cnt[u].find(v) == ct_cnt[u].end()) {
                ct_cnt[u].insert(make_pair(v,1));
            }else{
                ++ct_cnt[u][v];
            }
        }
    }

    // copy core time for old
    vector<vector<pair<int,int>>> core_t_old = core_t;

    // MAIN PROCESS
    int total_size = 0;
    int round_count = 0;
    vector<int> q(n_);
    int start = 0, end = 0;
    int t_s = ts_+1;
    int update_flag = false;
    int min_covered;
    #pragma omp parallel num_threads(threads_)
    {
        vector<vector<pair<int,int>>> mess(threads_);

        int tid = omp_get_thread_num();

        while (t_s <= te_) {
            #pragma omp master
            {
                min_covered = te_+1;
                start = 0;
                end = 0;
                // delete edges
                for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
                    int u = edges_[i].first;
                    int v = edges_[i].second;

                    if (invalid(u,core_t,init_tid) || invalid(v,core_t,init_tid)) continue;
                    if (!v_a[u]){   // process u
                        if (ct_cnt[u].find(v) != ct_cnt[u].end()) {
                            --ct_cnt[u][v];
                            if (ct_cnt[u][v] == 0) ct_cnt[u].erase(v);
                        }
                        if (ct_cnt[u].size()<k_){
                            q[end] = u;
                            end++;
                            v_a[u] = true;
                        }
                    }
                    if (!v_a[v]) {  // process v
                        if (ct_cnt[v].find(u) != ct_cnt[v].end()) {
                            --ct_cnt[v][u];
                            if (ct_cnt[v][u] == 0) ct_cnt[v].erase(u);
                        }
                        if (ct_cnt[v].size() < k_) {
                            q[end] = v;
                            end++;
                            v_a[v] = true;
                        }
                    }
                }
                if (end > start) update_flag = true;
            }
            #pragma omp barrier

            // process queue
            while (start < end){
                #pragma omp master
                {
                    total_size += end - start;
                    round_count += 1;
                }
                vector<bool> v_b(n_, false);
                vector<int> updated;

                #pragma omp for schedule(static)
                for (int j = start; j < end; j++) {     // process queue with parallel threads
                    int u = q[j];
                    updated.emplace_back(u);

                    // LocalCT
                    vector<int> nbr_t;
                    vector<int> bm_history;
                    int ct = 0;
                    for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                        int t = nbr_[u][i].second;
                        if (nbr_t.size() >= k_ && t > ct) break;
                        if (t < t_s){
                            offset[u] = i+1;
                            continue;
                        }

                        int v = nbr_[u][i].first;
                        if (invalid(v,core_t_old, init_tid) || v_b[v]) continue;
                        v_b[v] = true;
                        int v_t = core_t_old[v].back().second;
                        nbr_t.emplace_back(max(t,v_t));
                        bm_history.emplace_back(v);

                        if (nbr_t.size() <= k_) ct = max(ct,v_t);
                    }
                    for (auto &v:bm_history) v_b[v] = false;
                    int new_t = te_+1;
                    if (nbr_t.size() >= k_){
                        nth_element(nbr_t.begin(),nbr_t.begin()+k_-1,nbr_t.end());
                        new_t = nbr_t[k_-1];
                    }

                    // insert into index
                    int old_t = core_t[u].back().second;
                    if (old_t < min_covered) {
                        #pragma omp critical
                        {
                            min_covered = old_t;
                        }
                    }
                    if (core_t[u].back().first == t_s){
                        core_t[u].back().second = new_t;
                    }else{
                        core_t[u].emplace_back(make_pair(t_s,new_t));
                    }

                    // re-compute ct_cnt[u]
                    ct_cnt[u].clear();
                    vector<bool> v_c(n_, false);
                    for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                        int t = nbr_[u][i].second;
                        int v = nbr_[u][i].first;
                        if (t > new_t) break;
                        if (invalid(v,core_t_old,init_tid) || core_t_old[v].back().second > new_t) continue;

                        if (new_t != te_+1){
                            if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                                ct_cnt[u].insert(make_pair(v,1));
                            }else{
                                ++ct_cnt[u][v];
                            }
                        }

                        if (v_c[v]) continue;
                        v_c[v] = true;
                        if (new_t > core_t_old[v].back().second) {
                            mess[v/chunk_size_].emplace_back(make_pair(v,u));
                        }
                    }
                }

                // update ct_old and visied
                #pragma omp critical
                {   
                    for (int i = 0; i < mess.size(); ++i) {
                        for (auto vu : mess[i]) {
                            messages_[i].emplace_back(vu);
                        }
                    }
                    for (int i = 0; i < mess.size(); ++i) {
                        mess[i].clear();
                    }
                    for (int u : updated) {
                        core_t_old[u] = core_t[u];
                        v_a[u] = false;
                    }
                    start += updated.size();
                    updated.clear();
                }
                #pragma omp barrier

                for (auto e : messages_[omp_get_thread_num()]) {
                    int u = e.first;
                    int v = e.second;
            
                    if (v_a[u]) continue;

                    //   received > self
                    if (core_t[v].back().second > core_t[u].back().second) {    // if new CT[v] is later than CT[u] now, then v is not CT neighbor of u anymore
                        if (ct_cnt[u].find(v) != ct_cnt[u].end()) {
                            ct_cnt[u].erase(v);
                            if (ct_cnt[u].size() < k_) {
                                #pragma omp critical
                                {
                                    q[end] = u;
                                    end++;
                                    v_a[u] = true;  // can be not atomic, only this thread uses v_a[u]
                                }
                                // updated.emplace_back(u);
                            }
                        }
                    }
                }
                messages_[omp_get_thread_num()].clear();
                #pragma omp barrier
            }
            #pragma omp master
            {
                if (update_flag) {
                    t_cores(t_s, te_, min_covered, core_t);
                    update_flag = false;
                }
                t_s++;
            }
            #pragma omp barrier
        }
        #pragma omp barrier
    }
    cout << "average queue size = " << total_size / round_count << endl;

    #ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
        printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);

        struct rusage rUsage;
        getrusage(RUSAGE_SELF, &rUsage);
        long ms = rUsage.ru_maxrss;
        printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
    #endif

    // write_index(core_t);
}






// v4
void Graph::vcentric_parallel() {
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    // calculate chunk size
    if (t_ % threads_ == 0) {
        chunk_size_ = t_ / threads_;
    } else {
        chunk_size_ = (t_ + (threads_ - t_ % threads_)) / threads_;
    }

    // core decomposition
    int init_tid = 0;
    cores_.resize(threads_);

    // set start times
    vector<int> start_time;
    for (int i = 0; i * chunk_size_ < t_; ++i) {
        start_time.emplace_back(ts_ + i * chunk_size_);
    }

    // variables
    offsets_.resize(start_time.size());
    ct_cnts_.resize(start_time.size());
    core_ts_.resize(start_time.size());
    res_parallel_.resize(start_time.size());
    num_cores_parallel_.resize(start_time.size());

    for (int tid = 0; tid < start_time.size(); ++tid) {
        connect(tid, start_time[tid]);
        t_cores_parallel(start_time[tid], te_, -1, core_ts_[tid], tid);
    }

    // final result
    vector<vector<pair<int,int>>> core_t_final(n_);

    // MAIN PROCESS
    #pragma omp parallel for num_threads(threads_) schedule(static) ordered
    for (int stidx = 0; stidx < start_time.size(); stidx++) {
        // allocate start and end time for each thread
        int tid = omp_get_thread_num();
        int t_start = start_time[stidx];
        int t_end;
        if (stidx == start_time.size()-1) {
            t_end = te_+1;
        } else {
            t_end = start_time[stidx+1];
        }
        printf("thread ID = %d gets start = %d, end = %d\n", tid, t_start, t_end);

        // declare variables
        vector<bool> v_a(n_, false);
        vector<bool> v_b(n_, false);

        vector<int> offset = offsets_[stidx];
        vector<unordered_map<int,int>> ct_cnt = ct_cnts_[stidx];
        vector<vector<pair<int,int>>> core_t = core_ts_[stidx];

        // re-compute ct_cnt before start
        for (int u = 0; u < n_; ++u) {
            if (invalid(u, core_t, tid)) continue;
            ct_cnt[u].clear();

            int t = core_t[u].front().second;
            for (int i = offset[u]; i < nbr_[u].size(); ++i){
                if (nbr_[u][i].second > t) break;
                // if (i.second < t_start) continue;

                int v = nbr_[u][i].first;
                if (cores_[tid][v] < k_ || t < core_t[v].front().second) continue;
                if (ct_cnt[u].find(v) == ct_cnt[u].end()) {
                    ct_cnt[u].insert(make_pair(v,1));
                }else{
                    ++ct_cnt[u][v];
                }
            }
        }

        queue<int> q;
        bool update_flag = false;
        for (int t_s = t_start+1; t_s < t_end; ++t_s) {
            int min_covered = te_+1;
            for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
                int u = edges_[i].first;
                int v = edges_[i].second;

                if (invalid(u,core_t, tid) || invalid(v,core_t, tid)) continue;
                if (!v_a[u]){   // process u
                    if (ct_cnt[u].find(v) != ct_cnt[u].end()) {
                        --ct_cnt[u][v];
                        if (ct_cnt[u][v] == 0) ct_cnt[u].erase(v);
                    }
                    if (ct_cnt[u].size() < k_){
                        q.push(u);
                        v_a[u] = true;
                    }
                }
                if (!v_a[v]) {  // process v
                    if (ct_cnt[v].find(u) != ct_cnt[v].end()) {
                        --ct_cnt[v][u];
                        if (ct_cnt[v][u] == 0) ct_cnt[v].erase(u);
                    }
                    if (ct_cnt[v].size() < k_) {
                        q.push(v);
                        v_a[v] = true;
                    }
                }
            }
            if (!q.empty()) update_flag = true;

            // process queue
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                v_a[u] = false;

                // LocalCT
                vector<int> nbr_t;
                vector<int> bm_history;
                int ct = 0;
                for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                    int t = nbr_[u][i].second;
                    if (nbr_t.size() >= k_ && t > ct) break;
                    if (t < t_s){
                        offset[u] = i+1;
                        continue;
                    }

                    int v = nbr_[u][i].first;
                    if (invalid(v,core_t, tid) || v_b[v]) continue;
                    v_b[v] = true;
                    int v_t = core_t[v].back().second;
                    nbr_t.emplace_back(max(t,v_t));
                    bm_history.emplace_back(v);

                    if (nbr_t.size() <= k_) ct = max(ct,v_t);
                }
                for (auto &v:bm_history) v_b[v] = false;
                int new_t = te_+1;
                if (nbr_t.size() >= k_){
                    nth_element(nbr_t.begin(),nbr_t.begin()+k_-1,nbr_t.end());
                    new_t = nbr_t[k_-1];
                }

                // insert into index
                int old_t = core_t[u].back().second;
                min_covered = min(old_t, min_covered);
                if (core_t[u].back().first == t_s){
                    core_t[u].back().second = new_t;
                }else{
                    core_t[u].emplace_back(make_pair(t_s,new_t));
                }

                // re-compute ct_cnt[u]
                ct_cnt[u].clear();
                for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                    int t = nbr_[u][i].second;
                    int v = nbr_[u][i].first;
                    if (t > new_t) break;
                    if (invalid(v,core_t, tid) || core_t[v].back().second > new_t) continue;

                    if (new_t != te_+1){
                        if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                            ct_cnt[u].insert(make_pair(v,1));
                        }else{
                            ++ct_cnt[u][v];
                        }
                    }

                    // add neighbor to queue if necessary
                    if (v_a[v]) continue;
                    if (core_t[v].back().second < old_t || new_t <= core_t[v].back().second) continue;
                    ct_cnt[v].erase(u);
                    if (ct_cnt[v].size() < k_){
                        q.push(v);
                        v_a[v] = true;
                    }
                }
            }
            // compute distinct k-cores
            if (update_flag) {
                t_cores_parallel(t_s, te_, min_covered, core_t, stidx);
                update_flag = false;
            }
        }
        #pragma omp ordered
        {
            for (int u = 0; u < n_; u++) {
                for (int i = 0; i < core_t[u].size(); i++) {
                    if (!core_t_final[u].empty() && i == 0 && core_t_final[u].back().second == core_t[u][i].second) {
                        continue;
                    } else {
                        core_t_final[u].emplace_back(core_t[u][i]);
                    }
                }
            }
        }
    }

    #ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
        printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);

        struct rusage rUsage;
        getrusage(RUSAGE_SELF, &rUsage);
        long ms = rUsage.ru_maxrss;
        printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
    #endif

    // write_index(core_t_final);
    // write_res_parallel();
    // print_res_size();
    int total_num_cores = 0;
    for (auto v : num_cores_parallel_) total_num_cores += v;
    cout << "number of results = " << total_num_cores << endl;
}