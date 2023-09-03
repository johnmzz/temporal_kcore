#include "myGraph.h"
#include <omp.h>

Graph::Graph() {
    idx_size_ = 0;
    min_k_ = 2;
    nbr_cnt_ = nullptr;
    v_a_ = nullptr;
    v_b_ = nullptr;
    core_t_ = nullptr;
    t_offset_ = nullptr;
    core_ = nullptr;
    cd_ = nullptr;
    ct_cnt_ = nullptr;
    ctn_ = nullptr;
    ct_ = nullptr;
    // messages_ = nullptr;
}

Graph::~Graph() {
    delete[] nbr_cnt_;
    delete[] v_a_;
    delete[] v_b_;
    delete[] core_t_;
    delete[] core_;
    delete[] ctn_;
    delete[] ct_;
    // delete[] messages_;
    //delete[] cd_;   // same address as nbr_cnt_
    //delete[] ct_cnt_;  // same address as nbr_cnt_
}

void Graph::print_graph() {
    cout << "nbr_ = " << endl;
    for (int u = 0; u < nbr_.size(); u++) {
        cout << u << ": [";
        for (auto vt : nbr_[u]) {
            cout << "<" << vt.first << "," << vt.second << ">,";
        }
        cout << "]" << endl;
    }
}

void Graph::print_nbr_cnt_() {
    cout << "nbr_cnt_ = " << endl; 
    for (int i = 0; i < n_; i++) {
        cout << "u = " << i << ": [";
        for (auto j = nbr_cnt_[i].begin(); j != nbr_cnt_[i].end(); j++) {
            cout << "(" << j->first << "," << j->second << "),";
        }
        cout << "]" << endl;
    }
}

void Graph::print_nbr_time() {
    cout << "nbr_time_ = " << endl;
    for (int i = 0; i < n_; i++) {
        cout << i << ": [";
        for (auto t : nbr_time_[i]) {
            cout << t << ",";
        }
        cout << "]\n";
    }
    cout << "offset_ = " << endl;
    for (int i = 0; i < n_; i++) {
        cout << i << ": [";
        for (auto t : offset_[i]) {
            cout << t << ",";
        }
        cout << "]\n";
    }
}

void Graph::print_graph_size() {
    printf("Graph size: %.2f MB.\n",(float)edges_.size()*3*sizeof(int)/1024/1024);
}

void Graph::print_ct(vector<pair<int,int>>* core_t) {
    cout << "core_t = " << endl;
    for (int u = 0; u < n_; u++) {
        cout << "u = " << u << ": [";
        for (auto v : core_t[u]) {
            cout << "<" << v.first << "," << v.second << ">,";
        }
        cout << "]," << endl;
    }
}

void Graph::print_local_ct() {
    cout << "ct_ = [" << endl;
    for (int i = 0; i < n_; i++) {
        cout << std::setw(2) << i << ": [";        
        if (!ct_[i].empty()) {
            for (auto v : ct_[i]) {
                cout << std::setw(2) << "<" << v.first << "," << v.second << ">,";
            }
            // for (int j = 0; j < t_; j++) {
            //     cout << std::setw(2) << ct_[i][j] << ",";
            // }
        }
        cout << "]\n";
    }
}

void Graph::print_ctn(unordered_map<int, int>* ct_cnt) {
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

void Graph::print_queue(std::queue<int> q) {
  while (!q.empty())
  {
    std::cout << q.front() << " ";
    q.pop();
  }
  std::cout << std::endl;
}

void Graph::print_message() {
    cout << "messages = \n";
    for (int i = 0; i < threads_; i++) {
        cout << "thread id = " << i << ": [";
        for (auto v : messages_[i]) {
            cout << "<" << v.first << "," << v.second << ">,";
        }
        cout << endl;
    }
}

void Graph::write_index(vector<pair<int,int>>* core_t) {
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

void Graph::init_nbr_cnt() {
    if (nbr_cnt_!= nullptr) delete[] nbr_cnt_;

    effective_m_ = 0;
    max_effective_deg_ = 0;
    nbr_cnt_ = new unordered_map<int,int>[n_];

    for (int u = 0; u < n_; ++u) {
        for(auto &i : nbr_[u]) {
            if (nbr_cnt_[u].find(i.first) != nbr_cnt_[u].end()) {
                ++nbr_cnt_[u][i.first];
            } else {
                nbr_cnt_[u].insert(make_pair(i.first,1));
            }
        }
        if(nbr_cnt_[u].size() > max_effective_deg_) max_effective_deg_ = nbr_cnt_[u].size();
        effective_m_ += nbr_cnt_[u].size();
    }
    effective_m_ /= 2;
}

void Graph::init_nbr_time() {
    nbr_time_ = vector<vector<int>>(n_, vector<int>());
    offset_ = vector<vector<int>>(n_, vector<int>());
    for (int u = 0; u < n_; u++) {
        nbr_time_[u].emplace_back(-1);
    }

    for (int u = 0; u < n_; u++) {
        for (int i = 0; i < nbr_[u].size(); i++) {
            int t = nbr_[u][i].second;

            if (t != nbr_time_[u].back()) {
                nbr_time_[u].emplace_back(t);
                offset_[u].emplace_back(i);
            }
        }
    }
}

void Graph::load(const string &path) {
    printf("Graph path: %s\n", path.c_str());
    printf("Loading Graph\n");

    ifstream ifs(path);
    if (!ifs.is_open()) {
        cerr << "open file failed!" << endl;
        exit(1);
    }

    n_ = 0;
    m_ = 0;
    max_deg_ = 0;
    t_ = 0;
    t_min_ = LONG_MAX;

    int u,v;
    long ts;
    while (ifs.good() && !ifs.eof()) {
        char line[200];
        ifs.getline(line, 200);
        if (line[0] < '0' || line[0] > '9') continue;
        sscanf(line, "%d %d %ld", &u, &v, &ts);

        if (ts < t_min_) t_min_ = ts;

        if (u == v) continue;

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
        if (nbr_[u].size() > max_deg_) max_deg_ = nbr_[u].size();

        nbr_[v].emplace_back(make_pair(u,format_t));
        if (nbr_[v].size() > max_deg_) max_deg_ = nbr_[v].size();

        ++m_;
    }
    ifs.close();

    n_ = nbr_.size();
    t_ = t_new_to_old_.size();

    edges_idx_.emplace_back(edges_.size());

    init_nbr_cnt();

    // init_nbr_time();

    printf("n = %d, m = %d, effective_m = %d, max_deg = %d, max_effective_deg = %d.\n",n_,m_,effective_m_,max_deg_,max_effective_deg_);
    printf("span = %d, t_min = %ld\n",t_, t_min_);
    print_graph_size();

    if (v_a_ == nullptr) v_a_ = new bool[n_];
    if (v_b_ == nullptr) v_b_ = new bool[n_];
    for (int i = 0; i < n_; ++i) {
        v_a_[i] = false;
        v_b_[i] = false;
    }
}

void Graph::truncate(int ts, int te) {
    max_deg_ = 0;
    vector<vector<pair<int, int>>> nbr(n_);
    for (long t = ts; t <= te; t++) {
        for (int i = edges_idx_[t]; i < edges_idx_[t+1]; i++) {
            int u = edges_[i].first;
            int v = edges_[i].second;
            
            nbr[u].emplace_back(make_pair(v,t));
            if (nbr[u].size() > max_deg_) max_deg_ = nbr[u].size();

            nbr[v].emplace_back(make_pair(u,t));
            if (nbr[v].size() > max_deg_) max_deg_ = nbr[v].size();
        }
    }
    nbr_ = nbr;
    init_nbr_cnt();

    printf("After truncate:\n");
    printf("n = %d, m = %d, effective_m = %d, max_deg = %d, max_effective_deg = %d.\n",n_,m_,effective_m_,max_deg_,max_effective_deg_);
    printf("span = %d, t_min = %ld\n",t_, t_min_);
    print_graph_size();
}

void Graph::truncate_t(int ts, int te) {
    t_ = te - ts + 1;
    max_deg_ = 0;
    vector<vector<pair<int, int>>> nbr(n_);
    for (long t = ts; t <= te; t++) {
        for (int i = edges_idx_[t]; i < edges_idx_[t+1]; i++) {
            int u = edges_[i].first;
            int v = edges_[i].second;
            
            nbr[u].emplace_back(make_pair(v,t-ts));
            if (nbr[u].size() > max_deg_) max_deg_ = nbr[u].size();

            nbr[v].emplace_back(make_pair(u,t-ts));
            if (nbr[v].size() > max_deg_) max_deg_ = nbr[v].size();
        }
    }
    nbr_ = nbr;
    init_nbr_cnt();

    printf("After truncate:\n");
    printf("n = %d, m = %d, effective_m = %d, max_deg = %d, max_effective_deg = %d.\n",n_,m_,effective_m_,max_deg_,max_effective_deg_);
    printf("span = %d, t_min = %ld\n",t_, t_min_);
    print_graph_size();
}

void Graph::core_decomposition(int _k) {
    printf("starting core decomposition...\n");
    if (core_ == nullptr) core_ = new int[n_];

    t_offset_ = new int[n_];

    auto vert = new int[n_];
    auto bin = new int[max_effective_deg_+1];
    memset(bin,0,sizeof(int)*(max_effective_deg_+1));

    for (int u = 0; u < n_; ++u) {
        int d = nbr_cnt_[u].size();
        core_[u] = d;
        ++bin[d];
    }

    int offset = 0;
    for (int i = 0; i <= max_effective_deg_; ++i) {
        int num = bin[i];
        bin[i] = offset;
        offset += num;
    }

    for (int u = 0; u < n_; ++u) {
        t_offset_[u] = bin[core_[u]];
        vert[t_offset_[u]] = u;
        bin[core_[u]]++;
    }

    for (int i = max_effective_deg_; i >= 1; --i) bin[i] = bin[i - 1];
    bin[0] = 0;
    k_max_ = 0;

    int i;
    for (i = 0; i < n_; ++i) {
        int u = vert[i];
        if (core_[u] >= _k) {
            k_max_ = core_[u];
            break;
        }

        for (auto& item : nbr_[u]) {
            if (v_a_[item.first]) continue;
            v_a_[item.first] = true;
            if (core_[item.first] > core_[u]){
                int dv = core_[item.first], pv = t_offset_[item.first];
                int pw = bin[dv], w = vert[pw];
                if (item.first != w){
                    t_offset_[item.first] = pw, vert[pv] = w;
                    t_offset_[w] = pv, vert[pw] = item.first;
                }
                ++bin[dv];
                --core_[item.first];
            }
        }
        for (auto& item : nbr_[u]){
            v_a_[item.first] = false;
        }
        if (core_[u] > k_max_) k_max_ = core_[u];
    }
    for (;i < n_; i++) {
        int u = vert[i];
        core_[u] = _k;
    }
    delete[] bin;
    delete[] vert;
    delete[] t_offset_;
    cout << "Initial decomp complete.." << endl;
}

void Graph::compute_core_deg(const int &t_s) {
    if (cd_ == nullptr){
        cd_ = nbr_cnt_;
    }

    for (int u = 0; u < n_; ++u) {
        cd_[u].clear();
        for (int i = nbr_[u].size()-1;i>=0;--i){
            int v = nbr_[u][i].first;
            int t = nbr_[u][i].second;
            if (t < t_s) break;

            if (core_[v] < core_[u]) continue;

            if (cd_[u].find(v) == cd_[u].end()){
                cd_[u].insert(make_pair(v,1));
            }else{
                ++cd_[u][v];
            }
        }
    }
}

void Graph::init_core_time(int _ts, int _te, int _k) {
    int* core = new int[n_];
    for (int u = 0; u < n_; ++u) {
        core[u] = core_[u];
    }

    queue<int> q;
    int* cnt = new int[k_max_+1];
    for (int t_e = _te; t_e >= _ts; --t_e) {
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e+1]; ++i) { // delete edge
            int u = edges_[i].first;
            int v = edges_[i].second;

            // process u
            if (core[u] <= core[v]){
                --cd_[u][v];
                if (cd_[u][v] == 0){
                    cd_[u].erase(v);
                    if (cd_[u].size() < _k && !v_a_[u]){
                        q.push(u);
                        v_a_[u] = true;
                    }
                }
            }

            // process v
            if (core[v] <= core[u]){
                --cd_[v][u];
                if (cd_[v][u] == 0){
                    cd_[v].erase(u);
                    if (cd_[v].size() < _k && !v_a_[v]){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }
        }

        while (!q.empty()){
            int u = q.front();
            q.pop();
            int oc = core[u];
            memset(cnt,0,sizeof(int)*(oc+1));

            // LocalCore
            for (int i = 0; i < nbr_[u].size(); ++i) {
                int v = nbr_[u][i].first;
                int t = nbr_[u][i].second;
                if (t >= t_e) break;
                if (v_b_[v]) continue;
                v_b_[v] = true;

                ++cnt[core[v] < core[u] ? core[v]:core[u]];
            }
            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b_[nbr_[u][i].first] = false;
            }

            int cd = 0;
            for (int k = oc; k >= 0; --k) {
                cd += cnt[k];
                if(cd >= k){
                    core[u] = k;
                    break;
                }
            }

            // update cd;
            cd_[u].clear();
            for (int i = 0; i < nbr_[u].size(); ++i) {
                int v = nbr_[u][i].first;
                int t = nbr_[u][i].second;
                if (t >= t_e) break;
                if (core[v] < core[u]) continue;
                if (cd_[u].find(v) == cd_[u].end()) cd_[u].insert(make_pair(v,1));
                else ++cd_[u][v];
            }

            // add affected neighbor to the queue
            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                int v = nbr_[u][i].first;
                if (core[u] < core[v] && core[v] <= oc && !v_b_[v]){    // v affected if oldcore(u) > core(v) but now core(u) < core(v)
                    v_b_[v] = true;
                    cd_[v].erase(u);
                    if (!v_a_[v] && cd_[v].size() < _k){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }
            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b_[nbr_[u][i].first] = false;
            }

            if (core_[u] >= _k) {
                core_t_[u][_k].emplace_back(make_pair(0, t_e)); // if u popped from q, then its core changed at t_e, therefore put in Index
            }
        }
    }
    delete[] cnt;
    delete[] core;
}

void Graph::init_ct_cnt(int k) {
    // reuse nbr_cnt_ and avoid applying space
    ct_cnt_ = nbr_cnt_;

    for (int u = 0; u < n_; ++u) {
        if (invalid(u,k)) continue;
        t_offset_[u] = 0;
        ct_cnt_[u].clear();

        int t = core_t_[u][k].front().second;
        for (auto &i : nbr_[u]){
            if (i.second > t) break;
            int v = i.first;
            if (core_[v]<k || t < core_t_[v][k].front().second) continue;   // v's CT < u's CT, then v is in k-core earlier than u, thus add
            if (ct_cnt_[u].find(v) == ct_cnt_[u].end()){
                ct_cnt_[u].insert(make_pair(v,1));
            }else{
                ++ct_cnt_[u][v];
            }
        }
    }
}


void Graph::time_range_kcore(long _ts, long _te, int _k) {
    cout << "Query: ts = " << _ts << ", te = " << _te << ", k = " << _k << endl;

#ifdef _LINUX_
    struct timeval t_start, t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif

    int ts = t_old_to_new_[_ts];
    int te = t_old_to_new_[_te];
    cout << "after conversion: ts = " << ts << ", te = " << te << endl;

    truncate(ts, te);

    core_t_ = new vector<vector<pair<int,int>>>[n_];    // [u1: [k1:[pair1, pair2, ...], k2:[pair1, pair2, ...], ...], u2: [], ...]
    t_offset_ = new int[n_];
    
    printf("starting core decomposition...\n");
    core_decomposition(_k);
    for (int u = 0; u < n_; ++u) {
        core_t_[u].resize(core_[u]+1);      // for each u, set size of core_t_[u] to core(u)+1
    }
    printf("k_max = %d\n",k_max_);
    if (k_max_ < _k) {
        printf("queried k = %d exceed maximum core in G[ts,te] = %d\n", _k, k_max_);
        return;
    }

    printf("initialize core time.\n");
    compute_core_deg(ts);
    init_core_time(ts, te, _k);

    queue<int> q;
    for (int k = _k; k < _k+1; ++k) {
        cout << "k = " << k << endl;
        for (int i = 0; i < n_; i++) {
            v_a_[i] = false;
        }
        printf("Iteration k = %d.\n",k);
        init_ct_cnt(k);

        for (int t_s = ts+1; t_s <= te; ++t_s) {
            vector<int> cand;
            for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
                int u = edges_[i].first;
                int v = edges_[i].second;
                //cout << "delete edge = <" << u << "," << v << "," << t_s-1 << ">" << endl;

                if (invalid(u,k) || invalid(v,k)) continue;

                // process u
                if (!v_a_[u]){
                    del_nbr(u,v);
                    if (ct_cnt_[u].size()<k){
                        q.push(u);
                        v_a_[u] = true;
                    }
                }
                // process v
                if (!v_a_[v]) {
                    del_nbr(v, u);
                    if (ct_cnt_[v].size() < k) {
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }

            while (!q.empty()){
                int u = q.front();
                q.pop();
                v_a_[u] = false;

                ct_cnt_[u].clear();
                vector<int> nbr_t;
                vector<int> bm_history;
                int ct = 0;

                // LocalCT
                for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                    int t = nbr_[u][i].second;
                    if (nbr_t.size() >= k && t > ct) break;
                    if (t < t_s){
                        t_offset_[u] = i+1;
                        continue;
                    }

                    int v = nbr_[u][i].first;
                    if (invalid(v,k) || v_b_[v]) continue;
                    v_b_[v] = true; 
                    int v_t = core_t_[v][k].back().second;
                    nbr_t.emplace_back(max(t,v_t));
                    bm_history.emplace_back(v);

                    if (nbr_t.size() <= k) ct = max(ct,v_t);
                }
                for (auto &v:bm_history) v_b_[v] = false;

                int new_t = t_;
                if (nbr_t.size() >= k){
                    nth_element(nbr_t.begin(),nbr_t.begin()+k-1,nbr_t.end());
                    new_t = nbr_t[k-1];
                }

                // insert into index
                int old_t = core_t_[u][k].back().second;
                if (core_t_[u][k].back().first == t_s){
                    core_t_[u][k].back().second = new_t;
                }else{
                    core_t_[u][k].emplace_back(make_pair(t_s,new_t));
                }

                // re-compute ct_cnt_[u]
                for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                    int t = nbr_[u][i].second;
                    int v = nbr_[u][i].first;
                    if (t > new_t) break;
                    if (invalid(v,k) || core_t_[v][k].back().second > new_t) continue;

                    if (new_t != t_){
                        if (ct_cnt_[u].find(v) == ct_cnt_[u].end()){
                            ct_cnt_[u].insert(make_pair(v,1));
                        }else{
                            ++ct_cnt_[u][v];
                        }
                    }

                    // add neighbor to queue if necessary
                    if (v_a_[v]) continue;
                    if (core_t_[v][k].back().second < old_t || new_t <= core_t_[v][k].back().second) continue;
                    ct_cnt_[v].erase(u);
                    if (ct_cnt_[v].size() < k){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }
        }
        cout << "@@@@" << endl;
        exit(0);
    }
#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
    printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);


    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
#else
    clock_t end = clock();
    printf("Running time: %.2f s, %.2f min\n",(double)(end-start)/ CLOCKS_PER_SEC,(double)(end-start)/CLOCKS_PER_SEC/60);

#endif
    print_idx_size();
}




























void Graph::init_ct_v2(int ts, int te, int k, vector<int>& offset, vector<int>& ct_init) { // O(n*d)
    cout << "Initialize all CT for ts = " << ts << ", te = " << te << ", k = " << k << endl;

    vector<int> bm_history;
    vector<bool> visited(n_, false);

    for (int u = 1; u < n_; ++u) {
        if (core_[u] < k) continue;

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
            if (cnt == k) {
                ct_init[u] = t;
                break;
            }
        }
        for (auto &v : bm_history) visited[v] = false;
        bm_history.clear();
    }
    cout << "CT initialized." << endl;
}

void Graph::init_ctn_v2(int t_s, int t_e, int k, vector<int>& offset, vector<int>& ct_init, unordered_map<int, int>* ct_cnt) {

    for (int u = 1; u < n_; ++u) {
        if (core_[u] < k) continue;
        ct_cnt[u].clear();  // delete

        int t = ct_init[u];
        for (int i = offset[u]; i < nbr_[u].size(); ++i){
            // if (nbr_[u][i].second < t_s) continue;
            if (nbr_[u][i].second > t) break;

            int v = nbr_[u][i].first;
            if (core_[v] < k || t < ct_init[v]) continue;   // v's CT < u's CT, then v is in k-core earlier than u, thus add
            
            if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                ct_cnt[u].insert(make_pair(v,1));
            }else{
                ++ct_cnt[u][v];
            }
        }
    }
}

void Graph::local_ct_v2(int u, int t_s, int k, vector<bool> &visited, vector<int>& offset, vector<int>& ct_prev, vector<int>& ct_curr) {
    vector<int> nbr_t;
    vector<int> bm_history;
    int ub = 0;

    for (int i = offset[u]; i < nbr_[u].size(); ++i) {
        int t = nbr_[u][i].second;

        if (nbr_t.size() >= k && t > ub) break;

        // if (t < t_s) continue;

        int v = nbr_[u][i].first;
        if (core_[v] < k || visited[v] || ct_prev[v] == t_) continue;

        visited[v] = true;
        int v_t = ct_prev[v];
        int ct = max(t, v_t);

        nbr_t.emplace_back(ct);
        bm_history.emplace_back(v);

        if (nbr_t.size() <= k) ub = max(ub, ct);
    }
    if (nbr_t.size() >= k) {
        nth_element(nbr_t.begin(),nbr_t.begin()+k-1,nbr_t.end());
        ct_curr[u] = nbr_t[k-1];
    } else {
        ct_curr[u] = t_;
    }

    for (auto &v : bm_history) visited[v] = false;
}

void Graph::init_core_time_v2(int ts, int te, int k, int threads, vector<int>& offset, unordered_map<int, int>* ct_cnt, vector<pair<int,int>>* core_t) {
    printf("initialize core time.\n");
    vector<int> ct_old(n_, t_);

    init_ct_v2(ts, te, k, offset, ct_old);      // O(n*d)
    init_ctn_v2(ts, te, k, offset, ct_old, ct_cnt);     // O(n*d)

    vector<int> ct_new = ct_old;

    int round = 0;
    bool update = true;
    while (update) {
        update = false;
        round++;

        #pragma omp parallel num_threads(threads)
        {
            #pragma omp for schedule(static, chunk_size_)
            for (int u = 0; u < n_; ++u) {
                // cout << "thread " << omp_get_thread_num() << " processing u = " << u << endl;
                if (core_[u] < k || ct_old[u] == t_ || ct_cnt[u].size() >= k) continue;

                vector<bool> visited = vector<bool>(n_, false);
                
                int old_ct = ct_old[u];
                local_ct_v2(u, ts, k, visited, offset, ct_old, ct_new);
                int new_ct = ct_new[u];

                // re-compute ct_cnt[u]
                ct_cnt[u].clear();
                for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                    int t = nbr_[u][i].second;
                    // if (t < ts) continue;
                    if (t > new_ct) break;

                    int v = nbr_[u][i].first;
                    if (core_[v] < k) continue;

                    if (ct_old[v] > new_ct) continue;

                    if (new_ct != t_) {
                        if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                            ct_cnt[u].insert(make_pair(v,1));
                        }else{
                            ++ct_cnt[u][v];
                        }
                    }

                    // update neighbor
                    if (core_[v] < k || visited[v]) continue;
                    if (new_ct > ct_old[v]) {
                        visited[v] = true;
                        #pragma omp critical
                        {
                            messages_[v/chunk_size_].emplace_back(make_pair(v,u));
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

            #pragma omp for schedule(static, chunk_size_)
            for (int u = 0; u < n_; ++u) {
                if (ct_old[u] != ct_new[u]) {
                    ct_old[u] = ct_new[u];
                    update = true;
                }
            }
        }
    }

    for (int u = 1; u < n_; ++u) {
        if (core_[u] >= k) {
            core_t[u].emplace_back(make_pair(0, ct_new[u]));
        }
    }

    cout << "Initialization complete, round taken = " << round << endl;
}


void Graph::time_range_kcore_v2(long _ts, long _te, int k, int threads) {
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    // calculate chunk size
    threads_ = threads;
    if (n_ % threads == 0) {
        chunk_size_ = n_ / threads;
    } else {
        chunk_size_ = (n_ + (threads - n_ % threads)) / threads;
    }

    // transform time and truncate graph
    int ts = t_old_to_new_[_ts];
    int te = t_old_to_new_[_te];
    truncate(ts, te);
    
    // core decomposition
    core_decomposition(k);
    printf("k_max = %d\n",k_max_);
    if (k_max_ < k) {
        printf("queried k = %d exceed maximum core in G[ts,te] = %d\n", k, k_max_);
        return;
    }

    // declare variables
    vector<bool> v_a(n_, false);
    vector<int> offset(n_);
    unordered_map<int,int>* ct_cnt = new unordered_map<int,int>[n_];
    vector<pair<int,int>>* core_t = new vector<pair<int,int>>[n_];
    messages_ = vector<vector<pair<int, int>>>(threads);

    // initialize core time at ts = 0
    init_core_time_v2(0, te, k, threads, offset, ct_cnt, core_t);

    // copy core time for old
    vector<pair<int,int>>* core_t_old = new vector<pair<int,int>>[n_];
    for (int u = 0; u < n_; ++u) {
        core_t_old[u] = core_t[u];
    }

    // MAIN PROCESS
    auto chrono_start = chrono::high_resolution_clock::now();
    vector<int> q(n_);  // care: theoretically can exceed out of bound?
    for (int t_s = ts+1; t_s <= te; ++t_s) {
        int start = 0;
        int end = 0;

        // delete edges
        for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
            int u = edges_[i].first;
            int v = edges_[i].second;

            if (invalid_v2(u,k,core_t) || invalid_v2(v,k,core_t)) continue;
            if (!v_a[u]){   // process u
                if (ct_cnt[u].find(v) != ct_cnt[u].end()) {
                    --ct_cnt[u][v];
                    if (ct_cnt[u][v] == 0) ct_cnt[u].erase(v);
                }
                if (ct_cnt[u].size()<k){
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
                if (ct_cnt[v].size() < k) {
                    q[end] = v;
                    end++;
                    v_a[v] = true;
                }
            }
        }

        // process queue
        while (start < end){
            #pragma omp parallel num_threads(threads)   // create threads
            {   
                vector<bool> v_b(n_, false);
                vector<int> updated;

                #pragma omp for schedule(static)
                for (int j = start; j < end; j++) {     // process queue with parallel threads
                    int u = q[j];
                    #pragma omp critical
                    {
                        v_a[u] = false;
                    }

                    ct_cnt[u].clear();
                    vector<int> nbr_t;
                    vector<int> bm_history;
                    int ct = 0;

                    // LocalCT
                    for (int i = offset[u]; i < nbr_[u].size(); ++i) {  // TODO: local copy of offset
                        int t = nbr_[u][i].second;
                        if (nbr_t.size() >= k && t > ct) break;
                        if (t < t_s){
                            offset[u] = i+1;
                            continue;
                        }

                        int v = nbr_[u][i].first;
                        if (invalid_v2(v,k,core_t_old) || v_b[v]) continue;
                        v_b[v] = true;
                        int v_t = core_t_old[v].back().second;
                        nbr_t.emplace_back(max(t,v_t));
                        bm_history.emplace_back(v);

                        if (nbr_t.size() <= k) ct = max(ct,v_t);
                    }
                    for (auto &v:bm_history) v_b[v] = false;

                    int new_t = t_;
                    if (nbr_t.size() >= k){
                        nth_element(nbr_t.begin(),nbr_t.begin()+k-1,nbr_t.end());
                        new_t = nbr_t[k-1];
                    }

                    // insert into index
                    int old_t = core_t[u].back().second;
                    if (core_t[u].back().first == t_s){
                        core_t[u].back().second = new_t;
                    }else{
                        core_t[u].emplace_back(make_pair(t_s,new_t));
                    }
                    updated.emplace_back(u);    // no duplicates (each inner round)

                    // re-compute ct_cnt[u]
                    vector<bool> v_c(n_, false);
                    for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                        int t = nbr_[u][i].second;
                        int v = nbr_[u][i].first;
                        if (t > new_t) break;
                        if (invalid_v2(v,k,core_t_old) || core_t_old[v].back().second > new_t) continue;    //prove?

                        if (new_t != t_){
                            if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                                ct_cnt[u].insert(make_pair(v,1));
                            }else{
                                ++ct_cnt[u][v];
                            }
                        }

                        if (v_c[v]) continue;   // TODO: use v_b + bm_history?
                        v_c[v] = true;
                        
                        if (new_t > core_t_old[v].back().second) {
                            #pragma omp critical
                            {
                                messages_[v/chunk_size_].emplace_back(make_pair(v,u));
                            }
                        }
                    }
                }

                #pragma omp critical
                {
                    for (int u : updated) {
                        core_t_old[u] = core_t[u];  // TODO: only copy end?
                    }
                    start += updated.size();   // use updated.size()?
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
                            if (ct_cnt[u].size() < k) {
                                #pragma omp critical
                                {
                                    q[end] = u;
                                    end++;
                                    v_a[u] = true;  // can be not atomic, only this thread uses v_a[u]
                                }
                            }
                        }
                    }
                }
                messages_[omp_get_thread_num()].clear();    // TODO: use index
                #pragma omp barrier
            }
        }
    }
    auto chrono_stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(chrono_stop - chrono_start);
    cout << "duration = " << duration.count() << endl;

    #ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
        printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);

        struct rusage rUsage;
        getrusage(RUSAGE_SELF, &rUsage);
        long ms = rUsage.ru_maxrss;
        printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
    #endif

    write_index(core_t);

    delete[] ct_cnt;
    delete[] core_t;
    delete[] core_t_old;
}






























void Graph::time_range_kcore_v3(long _ts, long _te, int k, int threads) {
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    // calculate chunk size
    threads_ = threads;
    if (n_ % threads == 0) {
        chunk_size_ = n_ / threads;
    } else {
        chunk_size_ = (n_ + (threads - n_ % threads)) / threads;
    }

    // transform time and truncate graph
    int ts = t_old_to_new_[_ts];
    int te = t_old_to_new_[_te];
    truncate_t(ts, te);
    
    // core decomposition
    core_decomposition(k);
    printf("k_max = %d\n",k_max_);
    if (k_max_ < k) {
        printf("queried k = %d exceed maximum core in G[ts,te] = %d\n", k, k_max_);
        return;
    }

    // declare variables
    vector<bool> v_a(n_, false);
    vector<int> offset(n_);
    unordered_map<int,int>* ct_cnt = new unordered_map<int,int>[n_];
    vector<pair<int,int>>* core_t = new vector<pair<int,int>>[n_];
    messages_ = vector<vector<pair<int, int>>>(threads);

    // initialize core time at ts = 0
    init_core_time_v2(0, te, k, threads, offset, ct_cnt, core_t);

    // copy core time for old
    vector<pair<int,int>>* core_t_old = new vector<pair<int,int>>[n_];
    for (int u = 0; u < n_; ++u) {
        core_t_old[u] = core_t[u];
    }

    // MAIN PROCESS
    int total_size = 0;
    int round_count = 0;
    auto chrono_start = chrono::high_resolution_clock::now();
    vector<int> q(n_);  // care: theoretically can exceed out of bound?
    int start = 0, end = 0;
    int t_s = ts+1;
    #pragma omp parallel num_threads(threads)   // create threads
    {   
        auto chrono_start = chrono::high_resolution_clock::now();
        vector<vector<pair<int,int>>> mess(threads);
        while (t_s <= te) {
            #pragma omp master
            {   
                start = 0;
                end = 0;
                // delete edges
                for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
                    int u = edges_[i].first;
                    int v = edges_[i].second;

                    if (invalid_v2(u,k,core_t) || invalid_v2(v,k,core_t)) continue;
                    if (!v_a[u]){   // process u
                        if (ct_cnt[u].find(v) != ct_cnt[u].end()) {
                            --ct_cnt[u][v];
                            if (ct_cnt[u][v] == 0) ct_cnt[u].erase(v);
                        }
                        if (ct_cnt[u].size()<k){
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
                        if (ct_cnt[v].size() < k) {
                            q[end] = v;
                            end++;
                            v_a[v] = true;
                        }
                    }
                }
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
                    for (int i = offset[u]; i < nbr_[u].size(); ++i) {  // TODO: local copy of offset
                        int t = nbr_[u][i].second;
                        if (nbr_t.size() >= k && t > ct) break;
                        if (t < t_s){
                            offset[u] = i+1;
                            continue;
                        }

                        int v = nbr_[u][i].first;
                        if (invalid_v2(v,k,core_t_old) || v_b[v]) continue;
                        v_b[v] = true;
                        int v_t = core_t_old[v].back().second;
                        nbr_t.emplace_back(max(t,v_t));
                        bm_history.emplace_back(v);

                        if (nbr_t.size() <= k) ct = max(ct,v_t);
                    }
                    for (auto &v:bm_history) v_b[v] = false;
                    int new_t = t_;
                    if (nbr_t.size() >= k){
                        nth_element(nbr_t.begin(),nbr_t.begin()+k-1,nbr_t.end());
                        new_t = nbr_t[k-1];
                    }

                    // insert into index
                    int old_t = core_t[u].back().second;
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
                        if (invalid_v2(v,k,core_t_old) || core_t_old[v].back().second > new_t) continue;    //prove?

                        if (new_t != t_){
                            if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                                ct_cnt[u].insert(make_pair(v,1));
                            }else{
                                ++ct_cnt[u][v];
                            }
                        }

                        if (v_c[v]) continue;   // TODO: use v_b + bm_history?
                        v_c[v] = true;
                        if (new_t > core_t_old[v].back().second) {
                            // #pragma omp critical
                            // {
                            //     messages_[v/chunk_size_].emplace_back(make_pair(v,u));
                            // }
                            mess[v/chunk_size_].emplace_back(make_pair(v,u));
                        }
                    }
                }

                // updated ct_old and visited
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
                        core_t_old[u] = core_t[u];  // TODO: only copy end?
                        v_a[u] = false;
                    }
                    start += updated.size();   // use updated.size()?
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
                            if (ct_cnt[u].size() < k) {
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
                // #pragma omp barrier
                // #pragma omp critical 
                // {
                //     for (int u : updated) {
                //         if (!v_a[u]) {
                //             q[end] = u;
                //             end++;
                //             v_a[u] = true;
                //         }
                //     }
                // }
                messages_[omp_get_thread_num()].clear();    // TODO: use index
                #pragma omp barrier
            }
            #pragma omp master
            {
                t_s++;
            }
            #pragma omp barrier
        }
        #pragma omp barrier
        auto chrono_stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono_stop - chrono_start);
        #pragma omp critical
        {
            cout << "@@@, Thread = " << omp_get_thread_num() << " took " << duration.count() << endl;
        }
    }
    auto chrono_stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(chrono_stop - chrono_start);
    cout << "duration = " << duration.count() << endl;
    cout << "average queue size = " << total_size / round_count << endl;
    cout << "newest version" << endl;

    #ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
        printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);

        struct rusage rUsage;
        getrusage(RUSAGE_SELF, &rUsage);
        long ms = rUsage.ru_maxrss;
        printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
    #endif

    write_index(core_t);

    delete[] ct_cnt;
    delete[] core_t;
    delete[] core_t_old;
}































void Graph::init_ct_v3(int ts, int te, int k, vector<int>& offset, vector<int>& ct_init) { // O(n*d)
    cout << "Initialize all CT for ts = " << ts << ", te = " << te << ", k = " << k << endl;

    vector<int> bm_history;
    vector<bool> visited(n_, false);

    for (int u = 1; u < n_; ++u) {
        if (core_[u] < k) continue;

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
            if (cnt == k) {
                ct_init[u] = t;
                break;
            }
        }
        for (auto &v : bm_history) visited[v] = false;
        bm_history.clear();
    }
    cout << "CT initialized." << endl;
}

void Graph::init_ctn_v3(int t_s, int t_e, int k, vector<int>& offset, vector<int>& ct_init, unordered_map<int, int>* ct_cnt) {

    for (int u = 1; u < n_; ++u) {
        if (core_[u] < k) continue;
        ct_cnt[u].clear();  // delete

        int t = ct_init[u];
        for (int i = offset[u]; i < nbr_[u].size(); ++i){
            // if (nbr_[u][i].second < t_s) continue;
            if (nbr_[u][i].second > t) break;

            int v = nbr_[u][i].first;
            if (core_[v] < k || t < ct_init[v]) continue;   // v's CT < u's CT, then v is in k-core earlier than u, thus add
            
            if (ct_cnt[u].find(v) == ct_cnt[u].end()){
                ct_cnt[u].insert(make_pair(v,1));
            }else{
                ++ct_cnt[u][v];
            }
        }
    }
}


void Graph::local_ct_v3(int u, int t_s, int k, vector<bool> &visited, vector<int>& offset, vector<int>& ct_init) {
    vector<int> nbr_t;
    vector<int> bm_history;
    int ub = 0;

    for (int i = offset[u]; i < nbr_[u].size(); ++i) {
        int t = nbr_[u][i].second;

        if (nbr_t.size() >= k && t > ub) break;

        int v = nbr_[u][i].first;
        if (core_[v] < k || visited[v] || ct_init[v] == t_) continue;

        visited[v] = true;
        int v_t = ct_init[v];
        int ct = max(t, v_t);

        nbr_t.emplace_back(ct);
        bm_history.emplace_back(v);

        if (nbr_t.size() <= k) ub = max(ub, ct);
    }
    if (nbr_t.size() >= k) {
        nth_element(nbr_t.begin(),nbr_t.begin()+k-1,nbr_t.end());
        ct_init[u] = nbr_t[k-1];
    } else {
        ct_init[u] = t_;
    }

    for (auto &v : bm_history) visited[v] = false;
}



void Graph::init_core_time_v3(int ts, int te, int k, int threads, vector<int>& offset, unordered_map<int, int>* ct_cnt, vector<pair<int,int>>* core_t) {
    vector<int> ct_init(n_, t_);

    init_ct_v3(ts, te, k, offset, ct_init);      // O(n*d)
    init_ctn_v3(ts, te, k, offset, ct_init, ct_cnt);     // O(n*d)

    int round = 0;
    bool update = true;
    while (update) {
        update = false;

        for (int u = 1; u < n_; ++u) {
            if (core_[u] < k || ct_init[u] == t_ || ct_cnt[u].size() >= k) continue;
            
            vector<bool> visited(n_, false);

            int old_ct = ct_init[u];
            local_ct_v3(u, ts, k, visited, offset, ct_init);
            int new_ct = ct_init[u];

            // re-compte ct_cnt[u]
            for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                int t = nbr_[u][i].second;

                if (t > new_ct) break;

                int v = nbr_[u][i].first;
                if (core_[v] < k) continue;

                int v_ct = ct_init[v];
                if (v_ct > new_ct) continue;

                if (new_ct != t_) {
                    if (ct_cnt[u].find(v) == ct_cnt[u].end()) {
                        ct_cnt[u].insert(make_pair(v,1));
                    }else{
                        ++ct_cnt[u][v];
                    }
                }

                //update neighbor
                if (visited[v]) continue;
                if (v_ct < old_ct || new_ct <= v_ct) continue;

                ct_cnt[v].erase(u);
                visited[v] = true;
            }

            if (old_ct != new_ct) update = true;
        }
        round++;
    }

    for (int u = 1; u < n_; ++u) {
        if (core_[u] >= k) {
            core_t[u].emplace_back(make_pair(ts, ct_init[u]));
        }
    }

    cout << "Finished, round taken = " << round << endl;
}

void Graph::time_range_kcore_v4(long _ts, long _te, int k, int threads) {
    #ifdef _LINUX_
        struct timeval t_start, t_end;
        gettimeofday(&t_start, NULL);
    #endif

    // transform time and truncate graph
    int ts = t_old_to_new_[_ts];
    int te = t_old_to_new_[_te];
    truncate_t(ts, te);

    // calculate chunk size
    threads_ = threads;
    if (t_ % threads == 0) {
        chunk_size_ = t_ / threads;
    } else {
        chunk_size_ = (t_ + (threads - t_ % threads)) / threads;
    }
    
    // core decomposition
    core_decomposition(k);
    printf("k_max = %d\n",k_max_);
    if (k_max_ < k) {
        printf("queried k = %d exceed maximum core in G[ts,te] = %d\n", k, k_max_);
        return;
    }
    
    // set start times
    vector<int> start_time;
    for (int i = 0; i * chunk_size_ < t_; ++i) {
        start_time.emplace_back(i * chunk_size_);
    }

    // final result
    vector<pair<int,int>>* core_t_final = new vector<pair<int,int>>[n_];

    // MAIN PROCESS
    auto chrono_start = chrono::high_resolution_clock::now();
    #pragma omp parallel for num_threads(threads) schedule(static) ordered
    for (int stidx = 0; stidx < start_time.size(); stidx++) {
        // allocate start and end time for each thread
        int t_start = start_time[stidx];
        int t_end;
        if (stidx == start_time.size()-1) {
            t_end = t_;
        } else {
            t_end = start_time[stidx+1];
        }
        printf("thread ID = %d gets start = %d, end = %d\n", omp_get_thread_num(), t_start, t_end);

        // declare variables
        vector<bool> v_a(n_, false);
        vector<bool> v_b(n_, false);
        vector<int> offset(n_);
        unordered_map<int,int>* ct_cnt = new unordered_map<int,int>[n_];
        vector<pair<int,int>>* core_t = new vector<pair<int,int>>[n_];

        // run init
        init_core_time_v3(t_start, te, k, threads, offset, ct_cnt, core_t);

        // re-compute ct_cnt before start (necessary)
        for (int u = 0; u < n_; ++u) {
            if (invalid_v2(u,k,core_t)) continue;
            ct_cnt[u].clear();

            int t = core_t[u].front().second;
            for (int i = offset[u]; i < nbr_[u].size(); ++i){
                if (nbr_[u][i].second> t) break;
                // if (i.second < t_start) continue;

                int v = nbr_[u][i].first;
                if (core_[v]<k || t < core_t[v].front().second) continue;
                if (ct_cnt[u].find(v) == ct_cnt[u].end()) {
                    ct_cnt[u].insert(make_pair(v,1));
                }else{
                    ++ct_cnt[u][v];
                }
            }
        }
        
        // vector<int> q(n_);
        // int start = 0, end = 0;
        cout << "t_start = " << t_start << ", t_end = " << t_end << endl;
        queue<int> q;
        for (int t_s = t_start+1; t_s < t_end; ++t_s) {
            // start = 0;
            // end = 0;
            for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
                int u = edges_[i].first;
                int v = edges_[i].second;

                if (invalid_v2(u,k,core_t) || invalid_v2(v,k,core_t)) continue;
                if (!v_a[u]){   // process u
                    if (ct_cnt[u].find(v) != ct_cnt[u].end()) {
                        --ct_cnt[u][v];
                        if (ct_cnt[u][v] == 0) ct_cnt[u].erase(v);
                    }
                    if (ct_cnt[u].size()<k){
                        // q[end] = u;
                        // end++;
                        q.push(u);
                        v_a[u] = true;
                    }
                }
                if (!v_a[v]) {  // process v
                    if (ct_cnt[v].find(u) != ct_cnt[v].end()) {
                        --ct_cnt[v][u];
                        if (ct_cnt[v][u] == 0) ct_cnt[v].erase(u);
                    }
                    if (ct_cnt[v].size() < k) {
                        // q[end] = v;
                        // end++;
                        q.push(v);
                        v_a[v] = true;
                    }
                }
            }

            // process queue
            // while (start < end) {
            while (!q.empty()) {
                // int u = q[start];
                // start++;
                int u = q.front();
                q.pop();
                v_a[u] = false;

                // LocalCT
                vector<int> nbr_t;
                vector<int> bm_history;
                int ct = 0;
                for (int i = offset[u]; i < nbr_[u].size(); ++i) {
                    int t = nbr_[u][i].second;
                    if (nbr_t.size() >= k && t > ct) break;
                    if (t < t_s){
                        offset[u] = i+1;
                        continue;
                    }

                    int v = nbr_[u][i].first;
                    if (invalid_v2(v,k,core_t) || v_b[v]) continue;
                    v_b[v] = true;
                    int v_t = core_t[v].back().second;
                    nbr_t.emplace_back(max(t,v_t));
                    bm_history.emplace_back(v);

                    if (nbr_t.size() <= k) ct = max(ct,v_t);
                }
                for (auto &v:bm_history) v_b[v] = false;
                int new_t = t_;
                if (nbr_t.size() >= k){
                    nth_element(nbr_t.begin(),nbr_t.begin()+k-1,nbr_t.end());
                    new_t = nbr_t[k-1];
                }

                // insert into index
                int old_t = core_t[u].back().second;
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
                    if (invalid_v2(v,k,core_t) || core_t[v].back().second > new_t) continue;

                    if (new_t != t_){
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
                    if (ct_cnt[v].size() < k){
                        // q[end] = v;
                        // end++;
                        q.push(v);
                        v_a[v] = true;
                    }
                }
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

        delete[] ct_cnt;
        delete[] core_t;
    }
    auto chrono_stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(chrono_stop - chrono_start);
    cout << "duration = " << duration.count() << endl;

    #ifdef _LINUX_
        gettimeofday(&t_end, NULL);
        long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
        printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);

        struct rusage rUsage;
        getrusage(RUSAGE_SELF, &rUsage);
        long ms = rUsage.ru_maxrss;
        printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
    #endif

    write_index(core_t_final);
    delete[] core_t_final;
}