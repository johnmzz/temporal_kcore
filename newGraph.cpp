#include "newGraph.h"

void Graph::print_graph() {
    cout << "nbr_ = " << endl;
    for (int u = 0; u < nbr_.size(); u++) {
        if (!nbr_[u].empty()) {
            cout << u << ": [";
            for (auto neighbor : nbr_[u]) {
                cout << "<" << neighbor.v << "," << neighbor.t + 1 << "," << neighbor.eid + 1 << ">,";
            }            
            cout << "]" << endl;
        }
    }
}

void Graph::print_core() {
    cout << "[";
    for (int u = 0; u < core_.size(); u++) {
        cout << core_[u] << ",";
    }
    cout << "]\n";
}

void Graph::print_ct() {
    cout << "core_t_ = " << endl;
    for (int u = 0; u < n_; u++) {
        cout << "u = " << u << ": [";
        for (auto v : core_t_[u]) {
            cout << "<" << v.first+1 << "," << v.second+1 << ">,";
        }
        cout << "]," << endl;
    }
}

void Graph::print_ct_e() {
    cout << "ct_e_ = " << endl;
    for (int e = 0; e < m_; e++) {
        cout << "e = (" << edges_[e].first << "," << edges_[e].second << ") : [";
        for (auto i : ct_e_[e]) {
            cout << "<" << i.first+1 << "," << i.second+1 << ">,";
        }
        cout << "]," << endl;
    }
}

void Graph::print_skyline() {
    cout << "skyline_ = " << endl;
    for (int u = 0; u < skyline_.size(); u++) {
        
        if (version_ == "v")
            cout << "u = " << u << ": [";
        if (version_ == "e")
            cout << "(" << edges_[u].first << "," << edges_[u].second << ") : [";

        for (auto label : skyline_[u]) {
            cout << "<" <<label.v << ": " << label.start+1 << "," << label.end+1 << "," << label.act+1 << "," << label.last+1 << "," << label.aux << ">,";
        }
        cout << "]," << endl;
    }
}

void Graph::print_L() {
    cout << "L = [" << endl;
    for (auto label : L) {
        cout << "<" <<label.v << ": " 
             << label.start+1 << "," 
             << label.end+1 << "," 
             << label.act+1 << "," 
             << label.last+1 << "," 
             << label.aux << ","
             << label.prev << ","
             << label.next << ">,";
    }
    cout << "]\n";
}

void Graph::print_ba() {
    cout << "ba_ = " << endl;
    for (int t = 0; t < ba_.size(); t++) {
        cout << t+1 << ": [";
        for (auto i : ba_[t]) {
            cout << i << ",";
        }
        cout << "]\n";
    }
}

void Graph::print_bs() {
    cout << "bs_ = " << endl;
    for (int t = 0; t < bs_.size(); t++) {
        cout << t+1 << ": [";
        for (auto i : bs_[t]) {
            cout << i << ",";
        }
        cout << "]\n";
    }
}

void Graph::print_ll() {
    cout << "L = [";
    for (int h = 0; h != -1; h = L[h].next) {
        cout << h << ",";
    }
    cout << "]\n";
}

void Graph::output(int ts, int te, unordered_set<int>& R) {
    cout << "[" << ts+1 << "," << te+1 << "]: [";
    for (auto v : R) {
        if (version_ == "v")
            cout << v << ",";
        if (version_ == "e")
            cout << "(" << edges_[v].first << "," << edges_[v].second << "),";
    }
    cout << "]\n";

}

void Graph::ct_size() {
    ct_size_ = 0;
    cout << "core time size = ";
    for (int u = 0; u < n_; u++) {
        ct_size_ += core_t_[u].size();
    }
    cout << ct_size_ << endl;

    cout << "avg_deg = ";
    double total_deg = 0;
    for (int u = 0; u < n_; u++) {
        total_deg += nbr_[u].size();
    }
    double avg_deg = total_deg / (n_ - 1);
    cout << avg_deg << endl;

    cout << "ct_size * avg_deg = " << ct_size_ * avg_deg << endl;
}

void Graph::ct_e_size() {
    ct_size_ = 0;
    cout << "core time size = ";
    for (int e = 0; e < ct_e_.size(); e++) {
        ct_size_ += ct_e_[e].size();
    }
    cout << ct_size_ << endl;

    cout << "avg_deg = ";
    double total_deg = 0;
    for (int u = 0; u < n_; u++) {
        total_deg += nbr_[u].size();
    }
    double avg_deg = total_deg / (n_ - 1);
    cout << avg_deg << endl;

    cout << "ct_size * avg_deg = " << ct_size_ * avg_deg << endl;
}

string vectorToString(const vector<int>& vec) {
    if (vec.empty()) return "";

    ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i != 0) oss << ", ";
        oss << vec[i];
    }
    return oss.str();
}

void Graph::load(const string &path, long _ts, long _te, int k, string version, int write_res) {
    printf("Graph path: %s\n", path.c_str());
    write_res_ = write_res;
    version_ = version;

    ifstream ifs(path);
    if (!ifs.is_open()) {
        cerr << "open file failed!" << endl;
        exit(1);
    }

    n_ = 0;
    m_ = 0;
    k_ = k;
    num_res_ = 0;
    total_res_size_ = 0;
    
    int u,v;
    long ts;
    int edge_id;
    int format_t;
    
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
        edge_id = edges_.size() - 1;

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

        format_t = t_new_to_old_.size() - 1;

        // construct neighbor list
        Nbr n_u(v, format_t, edge_id);
        nbr_[u].emplace_back(n_u);

        Nbr n_v(u, format_t, edge_id);
        nbr_[v].emplace_back(n_v);
    }
    ifs.close();

    n_ = nbr_.size();
    m_ = edges_.size();
    t_ = format_t+1;

    edges_idx_.emplace_back(edges_.size());

    ts_ = t_old_to_new_[_ts];
    te_ = t_old_to_new_[_te];
    truncate(ts_, te_);
}

void Graph::truncate(int ts, int te) {
    // t_ = te - ts + 1;
    vector<vector<Nbr>> nbr(n_);
    for (long t = ts; t <= te; t++) {
        for (int i = edges_idx_[t]; i < edges_idx_[t+1]; i++) {
            int u = edges_[i].first;
            int v = edges_[i].second;

            Nbr n_u(v, t, i);
            nbr[u].emplace_back(n_u);

            Nbr n_v(u, t, i);
            nbr[v].emplace_back(n_v);
        }
    }
    nbr_ = nbr;
}

void Graph::core_decomposition(int ts) {
    int effective_m = 0;
    int max_effective_deg = 0;

    vector<unordered_map<int,int>> nbr_cnt(n_);

    for (int u = 0; u < n_; ++u) {
        for(auto &i : nbr_[u]) {
            if (i.t < ts) continue;
            if (nbr_cnt[u].find(i.v) != nbr_cnt[u].end()) {
                ++nbr_cnt[u][i.v];
            } else {
                nbr_cnt[u].insert(make_pair(i.v,1));
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
            if (v_a[item.v] || item.t < ts) continue;
            v_a[item.v] = true;
            if (core[item.v] > core[u]){
                int dv = core[item.v], pv = t_offset[item.v];
                int pw = bin[dv], w = vert[pw];
                if (item.v != w){
                    t_offset[item.v] = pw, vert[pv] = w;
                    t_offset[w] = pv, vert[pw] = item.v;
                }
                ++bin[dv];
                --core[item.v];
            }
        }

        for (auto& item : nbr_[u]){
            v_a[item.v] = false;
        }

        if (core[u] > k_max_) k_max_ = core[u];
    }
    cout << "k_max = " << k_max_ << endl;
    if (k_max_ < k_) {
        printf("queried k = %d exceed maximum core in G[ts,te] = %d\n", k_, k_max_);
        return;
    }

    core_ = core;
}

void Graph::compute_core_deg(const int &t_s, vector<unordered_map<int, int>>& cd, vector<int>& core) {
    for (int u = 0; u < n_; ++u) {
        cd[u].clear();
        for (int i = nbr_[u].size()-1;i>=0;--i){
            int v = nbr_[u][i].v;
            int t = nbr_[u][i].t;
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

void Graph::init_core_time(int _ts, int _te) {
    vector<int> core = core_;

    vector<unordered_map<int,int>> cd(n_);
    compute_core_deg(_ts, cd, core);

    vector<bool> v_a(n_, false);
    vector<bool> v_b(n_, false);
    vector<bool> v_e(m_, false);
    queue<int> q;
    vector<int> cnt(k_max_+1);
    for (int t_e = _te; t_e >= _ts; --t_e) {
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e+1]; ++i) { // delete edge
            int u = edges_[i].first;
            int v = edges_[i].second;

            if (!v_e[i] && core[u] >= k_ && core[v] >= k_) {
                ct_e_[i].emplace_back(make_pair(_ts, t_e));
                v_e[i] = true;
            }

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
                int v = nbr_[u][i].v;
                int t = nbr_[u][i].t;
                if (t >= t_e) break;
                if (t < _ts) continue;

                if (v_b[v]) continue;
                v_b[v] = true;

                ++cnt[core[v] < core[u] ? core[v]:core[u]];
            }
            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].t >= t_e) break;
                v_b[nbr_[u][i].v] = false;
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
                int v = nbr_[u][i].v;
                int t = nbr_[u][i].t;
                if (t >= t_e) break;
                if (t < _ts) continue;
                if (core[v] < core[u]) continue;
                if (cd[u].find(v) == cd[u].end()) cd[u].insert(make_pair(v,1));
                else ++cd[u][v];
            }

            // add affected neighbor to the queue
            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].t >= t_e) break;
                if (nbr_[u][i].t < _ts) continue;
                int v = nbr_[u][i].v;
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
                if (nbr_[u][i].t >= t_e) break;
                v_b[nbr_[u][i].v] = false;
            }

            if (oc >= k_ && core[u] < k_) {
                core_t_[u].emplace_back(make_pair(_ts, t_e)); // if u popped from q, then its core changed at t_e, therefore put in Index
                
                // compute core time of each edge in N(u)
                for (int i = 0; i < nbr_[u].size(); ++i) {
                    if (nbr_[u][i].t >= t_e) break;
                    if (nbr_[u][i].t < _ts) continue;

                    int v = nbr_[u][i].v;
                    int eid = nbr_[u][i].eid;

                    if (!v_e[eid] && core[v] >= k_) {
                        ct_e_[eid].emplace_back(make_pair(_ts, t_e));
                        v_e[eid] = true;
                    }
                }
            }
        }
    }
}

// compute core times for both vertices and edges
void Graph::e_core_time() {
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    core_.resize(n_);

    core_decomposition(ts_);
    if (k_max_ < k_) {
        printf("queried k = %d exceed maximum core in G[ts,te] = %d\n", k_, k_max_);
        return;
    }

    // declare variables
    vector<bool> v_a(n_, false);
    vector<bool> v_b(n_, false);
    vector<int> offset(n_);

    ct_cnt_.resize(n_);
    core_t_.resize(n_);

    ct_e_.resize(m_);

    // start-anchored core time
    init_core_time(ts_, te_);

    // compute ct_cnt before start
    for (int u = 0; u < n_; u++) {
        if (invalid(u)) continue;
        ct_cnt_[u].clear();

        int t = core_t_[u].front().second;   // the new core time of u
        for (int i = offset[u]; i < nbr_[u].size(); ++i){
            if (nbr_[u][i].t > t) break;
            if (nbr_[u][i].t < ts_) continue;

            int v = nbr_[u][i].v;
            if (core_[v] < k_ || t < core_t_[v].front().second) continue;
            if (ct_cnt_[u].find(v) == ct_cnt_[u].end()) {
                ct_cnt_[u].insert(make_pair(v,1));
            }else{
                ++ct_cnt_[u][v];
            }
        }
    }

    queue<int> q;
    for (int t_s = ts_+1; t_s <= te_; ++t_s) {
        for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
            int u = edges_[i].first;
            int v = edges_[i].second;

            if (ct_e_[i].empty()) continue;

            if (ct_e_[i].back().second != te_+1) {
                ct_e_[i].emplace_back(make_pair(t_s, te_+1));
            }

            if (invalid(u) || invalid(v)) continue;
            if (!v_a[u]){   // process u
                if (ct_cnt_[u].find(v) != ct_cnt_[u].end()) {
                    --ct_cnt_[u][v];
                    if (ct_cnt_[u][v] == 0) ct_cnt_[u].erase(v);
                }
                if (ct_cnt_[u].size() < k_){
                    q.push(u);
                    v_a[u] = true;
                }
            }
            if (!v_a[v]) {  // process v
                if (ct_cnt_[v].find(u) != ct_cnt_[v].end()) {
                    --ct_cnt_[v][u];
                    if (ct_cnt_[v][u] == 0) ct_cnt_[v].erase(u);
                }
                if (ct_cnt_[v].size() < k_) {
                    q.push(v);
                    v_a[v] = true;
                }
            }
        }

        while (!q.empty()){
            int u = q.front();
            q.pop();
            v_a[u] = false;

            // LocalCT
            vector<int> nbr_t;
            vector<int> bm_history;
            int ct = 0;
            for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                int t = nbr_[u][i].t;
                if (nbr_t.size() >= k_ && t > ct) break;
                if (t < t_s){   // ts_ must be less than t_s
                    offset[u] = i+1;
                    continue;
                }
               int v = nbr_[u][i].v;
                if (invalid(v) || v_b[v]) continue;
                v_b[v] = true;
                int v_t = core_t_[v].back().second;
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
            int old_t = core_t_[u].back().second;
            if (core_t_[u].back().first == t_s){
                core_t_[u].back().second = new_t;
            }else{
                core_t_[u].emplace_back(make_pair(t_s,new_t));
            }

            // re-compute ct_cnt[u]
            ct_cnt_[u].clear();
            for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                int t = nbr_[u][i].t;
                int v = nbr_[u][i].v;
                if (t > new_t) break;
                if (t < ts_) continue;  // should never be triggered as i starts at offset[u], which should be set to position > ts_
                if (invalid(v) || core_t_[v].back().second > new_t) continue;

                if (new_t != te_+1){     // if ct(u) invalid, no need to add into index, but need to update neighbors
                    if (ct_cnt_[u].find(v) == ct_cnt_[u].end()){
                        ct_cnt_[u].insert(make_pair(v,1));
                    }else{
                        ++ct_cnt_[u][v];
                    }
                }

                // add neighbor to queue if necessary
                if (v_a[v]) continue;
                if (core_t_[v].back().second < old_t || new_t <= core_t_[v].back().second) continue;
                ct_cnt_[v].erase(u);
                if (ct_cnt_[v].size() < k_){
                    q.push(v);
                    v_a[v] = true;
                }
            }

            for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                int t = nbr_[u][i].t;
                if (t > new_t) break;
                if (t < ts_) continue;

                int eid = nbr_[u][i].eid;
                if (ct_e_[eid].empty()) continue;

                if (new_t > ct_e_[eid].back().second) {
                    if (ct_e_[eid].back().first == t_s) {
                        ct_e_[eid].back().second = new_t;
                    }
                    else {
                        ct_e_[eid].emplace_back(make_pair(t_s, new_t));
                    }
                }
            }
        }
    }
    cout << "Finished compute CT." << endl;

    #ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
        printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);

        struct rusage rUsage;
        getrusage(RUSAGE_SELF, &rUsage);
        long ms = rUsage.ru_maxrss;
        printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
    #endif

    ct_e_size();
    // print_ct_e();
}

// only compute core times for vertices
void Graph::core_time() {
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    core_.resize(n_);
    ct_e_.resize(m_);

    core_decomposition(ts_);
    if (k_max_ < k_) {
        printf("queried k = %d exceed maximum core in G[ts,te] = %d\n", k_, k_max_);
        return;
    }

    // declare variables
    vector<bool> v_a(n_, false);
    vector<bool> v_b(n_, false);
    vector<int> offset(n_);

    ct_cnt_.resize(n_);
    core_t_.resize(n_);

    // start-anchored core time
    init_core_time(ts_, te_);

    // compute ct_cnt before start
    for (int u = 0; u < n_; u++) {
        if (invalid(u)) continue;
        ct_cnt_[u].clear();

        int t = core_t_[u].front().second;   // the new core time of u
        for (int i = offset[u]; i < nbr_[u].size(); ++i){
            if (nbr_[u][i].t > t) break;
            if (nbr_[u][i].t < ts_) continue;

            int v = nbr_[u][i].v;
            if (core_[v] < k_ || t < core_t_[v].front().second) continue;
            if (ct_cnt_[u].find(v) == ct_cnt_[u].end()) {
                ct_cnt_[u].insert(make_pair(v,1));
            }else{
                ++ct_cnt_[u][v];
            }
        }
    }

    queue<int> q;
    for (int t_s = ts_+1; t_s <= te_; ++t_s) {
        for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
            int u = edges_[i].first;
            int v = edges_[i].second;

            if (invalid(u) || invalid(v)) continue;
            if (!v_a[u]){   // process u
                if (ct_cnt_[u].find(v) != ct_cnt_[u].end()) {
                    --ct_cnt_[u][v];
                    if (ct_cnt_[u][v] == 0) ct_cnt_[u].erase(v);
                }
                if (ct_cnt_[u].size() < k_){
                    q.push(u);
                    v_a[u] = true;
                }
            }
            if (!v_a[v]) {  // process v
                if (ct_cnt_[v].find(u) != ct_cnt_[v].end()) {
                    --ct_cnt_[v][u];
                    if (ct_cnt_[v][u] == 0) ct_cnt_[v].erase(u);
                }
                if (ct_cnt_[v].size() < k_) {
                    q.push(v);
                    v_a[v] = true;
                }
            }
        }

        while (!q.empty()){
            int u = q.front();
            q.pop();
            v_a[u] = false;

            // LocalCT
            vector<int> nbr_t;
            vector<int> bm_history;
            int ct = 0;
            for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                int t = nbr_[u][i].t;
                if (nbr_t.size() >= k_ && t > ct) break;
                if (t < t_s){   // ts_ must be less than t_s
                    offset[u] = i+1;
                    continue;
                }
               int v = nbr_[u][i].v;
                if (invalid(v) || v_b[v]) continue;
                v_b[v] = true;
                int v_t = core_t_[v].back().second;
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
            int old_t = core_t_[u].back().second;
            if (core_t_[u].back().first == t_s){
                core_t_[u].back().second = new_t;
            }else{
                core_t_[u].emplace_back(make_pair(t_s,new_t));
            }

            // re-compute ct_cnt[u]
            ct_cnt_[u].clear();
            for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                int t = nbr_[u][i].t;
                int v = nbr_[u][i].v;
                if (t > new_t) break;
                if (t < ts_) continue;  // should never be triggered as i starts at offset[u], which should be set to position > ts_
                if (invalid(v) || core_t_[v].back().second > new_t) continue;

                if (new_t != te_+1){     // if ct(u) invalid, no need to add into index, but need to update neighbors
                    if (ct_cnt_[u].find(v) == ct_cnt_[u].end()){
                        ct_cnt_[u].insert(make_pair(v,1));
                    }else{
                        ++ct_cnt_[u][v];
                    }
                }

                // add neighbor to queue if necessary
                if (v_a[v]) continue;
                if (core_t_[v].back().second < old_t || new_t <= core_t_[v].back().second) continue;
                ct_cnt_[v].erase(u);
                if (ct_cnt_[v].size() < k_){
                    q.push(v);
                    v_a[v] = true;
                }
            }
        }
    }
    cout << "Finished compute CT." << endl;

    #ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
        printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);

        struct rusage rUsage;
        getrusage(RUSAGE_SELF, &rUsage);
        long ms = rUsage.ru_maxrss;
        printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
    #endif

    ct_size();
    // print_ct();
}








void Graph::baseline() {
    if (version_ == "v") {
        cout << "\n\n start compute all distinct vertex sets\n";
    }
    else if (version_ == "e") {
        cout << "\n\n start compute all distinct edge sets\n";
    }
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    ct_to_skyline();

    for (int ts = te_; ts >= ts_; ts--) {
        vector<vector<int>> bucket(t_);

        for (int u = 0; u < skyline_.size(); u++) {
            if (skyline_[u].empty()) continue;

            auto comp = [](const Label& a, int ts) { return a.start < ts; };

            auto it = lower_bound(skyline_[u].begin(), skyline_[u].end(), ts, comp);

            if (it != skyline_[u].end()) {
                bucket[it->end].emplace_back(u);
            }
        }

        vector<int> vset;
        for (int te = ts; te <= te_; te++) {
            if (bucket[te].empty()) continue;
            
            vset.insert(vset.end(), bucket[te].begin(), bucket[te].end());

            sort(vset.begin(), vset.end());
            string vset_str = vectorToString(vset);

            if (hres_.find(vset_str) == hres_.end()) {
                hres_[vset_str] = make_pair(ts, te);
                num_res_++;
                total_res_size_ += vset.size();

                if (write_res_) {
                    unordered_set<int> uset(vset.begin(), vset.end());
                    output(ts, te, uset); 
                }
            }
            else if (hres_[vset_str].second > te) {
                hres_[vset_str] = make_pair(ts, te);
                num_res_++;
                total_res_size_ += vset.size();

                if (write_res_) {
                    unordered_set<int> uset(vset.begin(), vset.end());
                    output(ts, te, uset);
                }
            }
        }
    }

    if (version_ == "v") {
        cout << "Finished compute vertex-sets, number of results = " << num_res_ << ", total res size = " << total_res_size_ << endl;
    }
    else if (version_ == "e") {
        cout << "Finished compute edge-sets, number of results = " << num_res_ << ", total res size = " << total_res_size_ << endl;
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
}












void Graph::ct_to_skyline() {
    if (version_ == "e") {        // use ct_e_
        skyline_.resize(m_);

        for (int u = 0; u < m_; u++) {
            for (int i = 0; i < ct_e_[u].size(); i++) {
                int ts, te;

                if (i == ct_e_[u].size()-1) {
                    if (ct_e_[u][i].second == te_+1) {
                        continue;
                    }
                    else {
                        ts = ct_e_[u][i].second;
                        te = ct_e_[u][i].second;
                        skyline_[u].emplace_back(Label(ts, te));
                    }
                }
                else {
                    ts = ct_e_[u][i+1].first-1;
                    te = ct_e_[u][i].second;
                    skyline_[u].emplace_back(Label(ts, te));
                }
            }
        }
    }
    else if (version_ == "v") {   // use core_t_
        skyline_.resize(n_);

        for (int u = 0; u < n_; u++) {
            for (int i = 0; i < core_t_[u].size(); i++) {
                int ts, te;

                if (i == core_t_[u].size()-1) {
                    if (core_t_[u][i].second == te_+1) {
                        continue;
                    }
                    else {
                        ts = core_t_[u][i].second;
                        te = core_t_[u][i].second;
                        skyline_[u].emplace_back(Label(ts, te));
                    }
                }
                else {
                    ts = core_t_[u][i+1].first-1;
                    te = core_t_[u][i].second;
                    skyline_[u].emplace_back(Label(ts, te));
                }
            }
        }
    }
    // print_skyline();
}

void Graph::init() {
    last_.resize(t_, -1);
    
    for (int u = 0; u < skyline_.size(); u++) {
        for (int i = 0; i < skyline_[u].size(); i++) {      // Labels ordered in start times
            // compute activation time
            if (i == 0) {
                skyline_[u][i].act = ts_;
            }
            else {
                skyline_[u][i].act = skyline_[u][i-1].start + 1;
            }

            // compute last valid time
            if (i == skyline_[u].size() - 1) {
                skyline_[u][i].last = te_;
            }
            else {
                skyline_[u][i].last = skyline_[u][i+1].end - 1;
            }

            skyline_[u][i].v = u;
            skyline_[u][i].aux = true;

            int t_last = skyline_[u][i].last;
            int t = skyline_[u][i].start;

            last_[t] = max(last_[t], t_last);
        }
    }
}

void Graph::compute_sets() {
    if (version_ == "v") {
        cout << "\n\n start compute all distinct vertex sets\n";
    }
    else if (version_ == "e") {
        cout << "\n\n start compute all distinct edge sets\n";
    }
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    ct_to_skyline();
    init();

    vector<vector<Label>> skyline_aux(skyline_.size());
    for (int u = 0; u < skyline_.size(); u++) {
        for (int i = 0; i < skyline_[u].size(); i++) {
            Label x = skyline_[u][i];
            x.end = x.last;
            x.aux = false;
            skyline_aux[u].emplace_back(x);
        }
    }
    
    // sort skyline labels
    vector<Label> labels;
    for (int u = 0; u < skyline_.size(); u++) {
        for (int i = 0; i < skyline_[u].size(); i++) {
            labels.emplace_back(skyline_[u][i]);
        }
    }
    std::sort(labels.begin(), labels.end(), [](const Label& a, const Label& b) {
        if (a.end == b.end && a.start == b.start) {
            return a.v < b.v;
        }
        else if (a.end == b.end) {
            return a.start < b.start;
        }
        return a.end < b.end;
    });

    // sort auxillary labels
    vector<Label> labels_aux;
    for (int u = 0; u < skyline_.size(); u++) {
        for (int i = 0; i < skyline_aux[u].size(); i++) {
            labels_aux.emplace_back(skyline_aux[u][i]);
        }
    }
    std::sort(labels_aux.begin(), labels_aux.end(), [](const Label& a, const Label& b) {
        if (a.end == b.end && a.start == b.start) {
            return a.v < b.v;
        }
        else if (a.end == b.end) {
            return a.start < b.start;
        }
        return a.end < b.end;
    });

    // merge labels and labels_aux ordered by end time, for the same end time labels always merged before labels_aux
    L.emplace_back(Label(-1,-1));   // insert dummy first
    int i = 0, j = 0;
    while (i < labels.size() && j < labels_aux.size()) {
        if (labels[i].end > labels_aux[j].end) {
            L.emplace_back(labels_aux[j]);
            j++;
        }
        else {
            L.emplace_back(labels[i]);
            i++;
        }
    }
    while (j < labels_aux.size()) {
        L.emplace_back(labels_aux[j]);
        j++;
    }

    // construct buckets
    ba_.resize(t_);
    bs_.resize(t_);
    for (int i = 1; i < L.size(); i++) {
        Label x = L[i];
        ba_[x.act].emplace_back(i);
        bs_[x.start].emplace_back(i);
    }

    // main process
    int head = 0;
    for (int t = ts_; t <= te_; t++) {
        if (t > ts_) {
            for (int x : bs_[t-1]) {
                remove(x);
            }
        }
        int h = head;

        for (int x : ba_[t]) {
            while (L[h].next != -1 && L[h].next < x) {
                h = L[h].next;
            }
            insert(x,h,L[h].next);
            h = x;
        }

        if (bs_[t].size() == 0) continue;

        // if (t == 2) print_L();

        enumerate(L[0].next, t, last_[t]);
    }

    if (version_ == "v") {
        cout << "Finished compute vertex-sets, number of results = " << num_res_ << ", total res size = " << total_res_size_ << endl;
    }
    else if (version_ == "e") {
        cout << "Finished compute edge-sets, number of results = " << num_res_ << ", total res size = " << total_res_size_ << endl;
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
}

void Graph::enumerate(int h, int t1, int t2) {
    int s = 0;
    unordered_set<int> R;

    while (h != -1) {
        if (L[h].end > t2) break;

        if (L[h].aux == false) {
            if (L[h].start == t1) s--;
            h = L[h].next;
            continue;
        }

        if (L[h].start == t1) s++;

        R.insert(L[h].v);
            
        if (s == 0) {
            h = L[h].next;
            continue;
        }

        int next = L[h].next;
        if (next != -1 && L[h].end == L[next].end && L[next].aux == true) {
            h = next;
            continue;
        }
        
        // collect result
        num_res_++;
        total_res_size_ += R.size();
        if (write_res_) output(t1, L[h].end, R);

        h = next;
    }
}

void Graph::remove(int x) {
    Label l = L[x];
    L[l.prev].next = l.next;
    if (l.next != -1) L[l.next].prev = l.prev;
}

void Graph::insert(int x, int a, int b) {
    L[x].next = b;
    L[x].prev = a;
    L[a].next = x;
    if (b != -1) L[b].prev = x;
}