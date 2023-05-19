#include "Graph.h"


Graph::Graph() {

    nbr_cnt_ = nullptr;
    ct_cnt_ = nullptr;
    min_k_ = 2;
    idx_size_ = 0;
    query_v_ = nullptr;
    query_deg_ = nullptr;
    core_t_ = nullptr;
    v_a_ = nullptr;
    v_b_ = nullptr;
    core_ = nullptr;
    t_offset_ = nullptr;
    cd_ = nullptr;
    log_f_ = nullptr;
}

Graph::~Graph() {

    delete[] nbr_cnt_;
    delete[] core_t_;

    delete[] core_;

//    delete[] ct_cnt_;
    delete[] t_offset_;

    delete[] v_a_;
    delete[] v_b_;
//    delete[] cd_;

    delete[] query_v_;
    delete[] query_deg_;
    if(log_f_!= nullptr) fclose(log_f_);
}

void Graph::load(const string &path, const int &total_edge_num, bool timestamp_third) {
    if(log_f_ != nullptr) fprintf(log_f_,"Graph path: %s, total edge count: %d\n",path.c_str(),total_edge_num);
    printf("Graph path: %s\n",path.c_str());
    printf("Loading Graph\n");

    ifstream ifs(path);
    if(!ifs.is_open()){
        cerr << "open file failed!" << endl;
        exit(-1);
    }


    n_ = 0;
    m_ = 0;
    max_deg_ = 0;
    t_ = 0;


    int u,v;
    long ts;
    while(ifs.good() && !ifs.eof()){
        char line[200];
        ifs.getline(line,200);
        if( line[0] < '0' || line[0] > '9' ) continue;
        if(timestamp_third){
            sscanf(line,"%d %d %ld",&u, &v, &ts);
        }else{
            int weight;
            sscanf(line,"%d %d %d %ld",&u, &v, &weight, &ts);
        }

        if(u == v) continue;

        edges_.emplace_back(make_pair(u,v));

        //adjust size of neighbor list if necessary.
        if(u+1 > nbr_.size()){
            nbr_.resize(u+1);
//            reverse_nbr_idx_.resize(u + 1);
        }
        if(v+1 > nbr_.size()){
            nbr_.resize(v + 1);
//            reverse_nbr_idx_.resize(v + 1);
        }

        long pre_ts = -1;
        if (!t_new_to_old_.empty()){
            pre_ts = t_new_to_old_.back();
        }

        if (ts != pre_ts){
            t_new_to_old_.emplace_back(ts);
            t_old_to_new_.insert(make_pair(ts,t_new_to_old_.size()-1));
            edges_idx_.emplace_back(edges_.size() - 1);
        }

        int format_t = t_new_to_old_.size()-1;

//        construct neighbor list
        nbr_[u].emplace_back(make_pair(v,format_t));
        if(nbr_[u].size() > max_deg_) max_deg_ = nbr_[u].size();

        nbr_[v].emplace_back(make_pair(u,format_t));
        if(nbr_[v].size() > max_deg_) max_deg_ = nbr_[v].size();

//        construct reverse neighbor index
//        reverse_nbr_idx_[u].emplace_back(nbr_[v].size() - 1);
//        reverse_nbr_idx_[v].emplace_back(nbr_[u].size() - 1);

        ++m_;
        if (total_edge_num!=-1 && m_ >= total_edge_num) break;
    }
    ifs.close();

    n_ = nbr_.size();
    t_ = t_new_to_old_.size();

    edges_idx_.emplace_back(edges_.size());

//    cout<<"Vertex Num: "<<n_<<", Edge Num: "<<m_<<"; Time Span: "<<t_<<endl;

    init_nbr_cnt();

    printf("n = %d, m = %d, effective_m = %d, max_deg = %d, max_effective_deg = %d.\n",n_,m_,effective_m_,max_deg_,max_effective_deg_);
    printf("span = %d.\n",t_);

    if(log_f_!= nullptr){
        fprintf(log_f_, "n = %d, m = %d, effective_m = %d, max_deg = %d, max_effective_deg = %d.\n",n_,m_,effective_m_,max_deg_,max_effective_deg_);
        fprintf(log_f_, "span = %d.\n",t_);
    }

    print_graph_size();


    if (v_a_ == nullptr) v_a_ = new bool[n_];
    if (v_b_ == nullptr) v_b_ = new bool[n_];
    for (int i = 0; i < n_; ++i) {
        v_a_[i] = false;
        v_b_[i] = false;
    }
    query_v_ = new vector<int>[n_];
    query_deg_ = new int[n_];

}

//initialize effective_deg_[], forward_nbr_cnt_[], and backward_nbr_cnt_[]
void Graph::init_nbr_cnt() {

    effective_m_ = 0;

    nbr_cnt_ = new unordered_map<int,int>[n_];


    max_effective_deg_ = 0;

    for (int u = 0; u < n_; ++u) {

        for(auto &i : nbr_[u]){
            if (nbr_cnt_[u].find(i.first) != nbr_cnt_[u].end()){
                ++nbr_cnt_[u][i.first];
            } else{
                nbr_cnt_[u].insert(make_pair(i.first,1));
            }
        }

        if(nbr_cnt_[u].size() > max_effective_deg_) max_effective_deg_ = nbr_cnt_[u].size();
        effective_m_ += nbr_cnt_[u].size();
    }

    effective_m_ /= 2;
}

void Graph::index() {

#ifdef _LINUX_
    struct timeval t_start,t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif


    core_t_ = new vector<vector<pair<int,int>>>[n_];    // [u1: [k1:[pair1, pair2, ...], k2:[pair1, pair2, ...], ...], u2: [], ...]
    t_offset_ = new int[n_];
    if (v_a_ == nullptr) v_a_ = new bool[n_];
    if (v_b_ == nullptr) v_b_ = new bool[n_];


    printf("starting core decomposition...\n");
//    compute core number for all edges
    core_decomposition();
    for (int u = 0; u < n_; ++u) {
        v_a_[u] = false;
        v_b_[u] = false;
        core_t_[u].resize(core_[u]+1);      // for each u, set size of core_t_[u] to core(u)+1
    }
    printf("k_max = %d\n",k_max_);

    printf("initialize core time.\n");
    compute_core_deg(0);
    init_core_time();


//    init ct_deg_ and ct_cnt
    queue<int> q;
    for (int k = min_k_; k <= k_max_; ++k) {
        printf("Iteration k = %d.\n",k);
        init_ct_cnt(k);

        for (int t_s = 1; t_s < t_; ++t_s) {
            vector<int> cand;
            for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
                int u = edges_[i].first;
                int v = edges_[i].second;
                // cout << "delete edge = <" << u << "," << v << "," << t_s-1 << ">" << endl;

                if (invalid(u,k) || invalid(v,k)) continue;

//                process u
                if (!v_a_[u]){
                    del_nbr(u,v);
                    if (ct_cnt_[u].size()<k){
                        q.push(u);
                        v_a_[u] = true;
                    }
                }


//                process v
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

//                  compute new core time of u
                for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {

                    int t = nbr_[u][i].second;
                    if (nbr_t.size() >= k && t > ct) break;
                    if (t < t_s){
                        t_offset_[u] = i+1;
                        continue;
                    }

                    int v = nbr_[u][i].first;
//                    v is not valid or has already been visited
                    if (invalid(v,k) || v_b_[v]) continue;

//                    mark v has been visited
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
//                    sort(nbr_t.begin(),nbr_t.end());
                    new_t = nbr_t[k-1];
                }

                int old_t = core_t_[u][k].back().second;
                if (core_t_[u][k].back().first == t_s){
                    core_t_[u][k].back().second = new_t;
                }else{
                    core_t_[u][k].emplace_back(make_pair(t_s,new_t));
                }

//                compute ct_cnt_[u] and add neighbor to queue if necessary
                for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
//                    compute ct_cnt_[u]
                    int t = nbr_[u][i].second;
                    int v = nbr_[u][i].first;
                    if (t > new_t) break;
                    if (invalid(v,k) || core_t_[v][k].back().second > new_t) continue;

                    if (new_t != t_){
                        if (ct_cnt_[u].find(v)==ct_cnt_[u].end() ){
                            ct_cnt_[u].insert(make_pair(v,1));
                        }else{
                            ++ct_cnt_[u][v];
                        }
                    }

//                    add neighbor to queue if necessary
                    if (v_a_[v]) continue;
                    if (core_t_[v][k].back().second < old_t || new_t <= core_t_[v][k].back().second) continue;
//                    del_nbr(v,u);
                    ct_cnt_[v].erase(u);
                    if (ct_cnt_[v].size() < k){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }

        }
    }





#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
    printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);
    if(log_f_ != nullptr) fprintf(log_f_,"Indexing time: %lld s\n",t_msec/1000);


    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    printf("Memory usage = %ldKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024);
    if(log_f_ != nullptr) fprintf(log_f_,"Memory usage = %ldKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024);
#else
    clock_t end = clock();
    printf("Running time: %.2f s, %.2f min\n",(double)(end-start)/ CLOCKS_PER_SEC,(double)(end-start)/CLOCKS_PER_SEC/60);

#endif
    if(log_f_ != nullptr) fprintf(log_f_,"kmax = %d\n",k_max_);
    print_idx_size();

    //print_index();
}



//core decomposition for all edges
void Graph::core_decomposition() {
    if (core_ == nullptr) core_ = new int[n_];

    auto vert = new int[n_];
    auto bin = new int[max_effective_deg_+1];  
    memset(bin,0,sizeof(int)*(max_effective_deg_+1));

    for (int u = 0; u < n_; ++u) {
        int d = nbr_cnt_[u].size();
        core_[u] = d;
        ++bin[d];
    }

    // sort based on degree
    int offset = 0;
    for (int i = 0; i <= max_effective_deg_ ; ++i) {
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

    for (int i = 0; i < n_; ++i) {
        int u = vert[i];

        for (auto& item : nbr_[u]){
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
        // early termination
        // if (core_[u] > _k) {
        //     core_[u] = k;
        //     for each core_[u] remaining in vert:
        //         core_u[u] = k;
        //     break;
        // }
    }

    delete[] bin;
    delete[] vert;
}

void Graph::init_core_time() {
//    vector<int> core_history;
    int* core = new int[n_];
    for (int u = 0; u < n_; ++u) {
        core[u] = core_[u];
    }

    queue<int> q;
    int* cnt = new int[k_max_+1];
    for (int t_e = t_-1; t_e >= 0 ; --t_e) {
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e + 1]; ++i) {

            int u = edges_[i].first;
            int v = edges_[i].second;
            //cout << "delete edge = <" << u << "," << v << "," << t_e << ">" << endl;

    //        process u (delete edge)
            if (core[u] <= core[v]){
                --cd_[u][v];
                if (cd_[u][v] == 0){
                    cd_[u].erase(v);
                    if (cd_[u].size() < core[u] && !v_a_[u]){
                        q.push(u);
                        v_a_[u] = true;
                    }
                }
            }

//            process v (delete edge)
            if (core[v] <= core[u]){
                --cd_[v][u];
                if (cd_[v][u] == 0){
                    cd_[v].erase(u);
                    if (cd_[v].size() < core[v] && !v_a_[v]){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }
        }



        while (!q.empty()){
            int u = q.front();
            q.pop();
            v_a_[u] = false;

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
            for (int k = oc; k >= 0 ; --k) {
                cd += cnt[k];
                if(cd >= k){
                    core[u] = k;
                    break;
                }
            }

//            update cd_; (redundant?)
            cd_[u].clear();
            for (int i = 0; i < nbr_[u].size(); ++i) {
                int v = nbr_[u][i].first;
                int t = nbr_[u][i].second;
                if (t >= t_e) break;
                if (core[v] < core[u]) continue;
                if (cd_[u].find(v) == cd_[u].end()) cd_[u].insert(make_pair(v,1));
                else ++cd_[u][v];
            }


    //        add influenced neighbor to the queue
            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                int v = nbr_[u][i].first;
                if (core[u] < core[v] && core[v] <= oc && !v_b_[v]){
                    v_b_[v] = true;
                    cd_[v].erase(u);
                    if (!v_a_[v] && cd_[v].size() < core[v]){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }

            for (int i = 0; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b_[nbr_[u][i].first] = false;
            }

            for (int k = oc; k > core[u]; --k) {
                core_t_[u][k].emplace_back(make_pair(0, t_e));  // if u popped from q, then its core changed at t_e, therefore put in Index
            }

        }

    }



    delete[] cnt;
    delete[] core;
}

void Graph::init_ct_cnt(int k) {

//    reuse nbr_cnt_ and avoid applying space
    ct_cnt_ = nbr_cnt_;

    for (int u = 0; u < n_; ++u) {
        if (invalid(u,k)) continue;
        t_offset_[u] = 0;
        ct_cnt_[u].clear();

        int t = core_t_[u][k].front().second;
        for (auto &i : nbr_[u]){
            if (i.second > t) break;
            int v = i.first;
            if (core_[v]<k || t < core_t_[v][k].front().second) continue;
            if (ct_cnt_[u].find(v) == ct_cnt_[u].end()){
                ct_cnt_[u].insert(make_pair(v,1));
            }else{
                ++ ct_cnt_[u][v];
            }
        }
    }
}

void Graph::test() {

    for (int u = 0; u < 20; ++u) {
        printf("vertex %d:\n",u);

        for (int k = min_k_; k < core_t_[u].size(); ++k) {
            printf("k=%d:\n",k);
            for (auto &i:core_t_[u][k]) printf("[%d,%d], ",i.first,i.second);
            printf("\n");
        }
        printf("\n");
    }

//    auto fp = fopen(R"(C:\Users\DW\Desktop\idx\email-b)","rb");
//    int c,d;
//    fread(&c,sizeof(unsigned int),1,fp);
//
//
//    for (int u = 0; u < n_; ++u) {
//        int cs;
//        fread(&cs,sizeof(int),1,fp);
//        if (cs != core_t_[u].size()) printf("neq %d\n",u);
//
//        for (int k = min_k_; k < cs; ++k) {
//            int ccs;
//            fread(&ccs,sizeof(int),1,fp);
//            if (ccs != core_t_[u][k].size()) printf("neq ccs %d[%d]\n",u,k);
//            for (int i = 0; i < ccs; ++i) {
//                int a,b;
//                fread(&a,sizeof(int),1,fp);
//                fread(&b,sizeof(int),1,fp);
//
//            }
//        }
//    }
//    fclose(fp);

}

void Graph::print_idx_size() {
    if (idx_size_ != 0) printf("Index size: %lld MB.",idx_size_/1024/1024);

    idx_size_ += sizeof(int);
    idx_size_ += sizeof(int)*n_;

    double average_t = 0;
    int average_t_d = 0;
    int max_t = 0;

    for (int u = 0; u < n_; ++u) {
        if (core_t_[u].size() <= min_k_) continue;
        for (int k = min_k_; k < core_t_[u].size(); ++k){
            idx_size_ += core_t_[u][k].size()*2*sizeof(int);
            average_t += core_t_[u][k].size();
            if (core_t_[u][k].size() > max_t) max_t = core_t_[u][k].size();
        }
        average_t_d += core_t_[u].size()-min_k_;
    }

    printf("Index size: %.2f MB.\n",(float)idx_size_/1024/1024);
    printf("Average T = %.2f, max T = %d.\n",average_t/average_t_d,max_t);

    if(log_f_ != nullptr) fprintf(log_f_,"Index size: %.2f MB\n",(float)idx_size_/1024/1024);
    if(log_f_ != nullptr) fprintf(log_f_,"Average T = %.2f, max T = %d.\n",average_t/average_t_d,max_t);
}

void Graph::print_graph_size() {
    printf("Graph size: %.2f MB.\n",(float)edges_.size()*3*sizeof(int)/1024/1024);
}

bool Graph::query(int u, int t_s, int t_e, int k) {
    if (core_t_[u].size() < k+1) return false;
    if (k < 2){
        printf("Enter a parameter larger than 1\n");
        return false;
    }
    auto it = upper_bound(core_t_[u][k].begin(),core_t_[u][k].end(),make_pair(t_s,t_e),cmp);
    --it;

    return it->second <= t_e;
}

void Graph::write_idx(const string &path) {
    auto fp = fopen(path.c_str(),"wb");
    fwrite(&n_,sizeof(unsigned int),1,fp);

    for (int u = 0; u < n_; ++u) {
        int cs = core_t_[u].size();
        fwrite(&cs,sizeof(int),1,fp);
        for (int k = min_k_; k < cs; ++k) {
            int ccs = core_t_[u][k].size();
            fwrite(&ccs,sizeof(int),1,fp);
            for (auto &i:core_t_[u][k]){
                fwrite(&i.first,sizeof(int),1,fp);
                fwrite(&i.second,sizeof(int),1,fp);
            }
        }
    }
    fclose(fp);
}

void Graph::load_idx(const string &path) {
    auto fp = fopen(path.c_str(),"rb");
    fread(&n_,sizeof(unsigned int),1,fp);
    k_max_ = 0;
    t_ = 0;

    if (core_t_ == nullptr) core_t_ = new vector<vector<pair<int,int>>>[n_];

    for (int u = 0; u < n_; ++u) {
        int cs;
        fread(&cs,sizeof(int),1,fp);
        core_t_[u].resize(cs);
        if(k_max_ < cs-1) k_max_ = cs-1;

        for (int k = min_k_; k < cs; ++k) {
            int ccs;
            fread(&ccs,sizeof(int),1,fp);

            for (int i = 0; i < ccs; ++i) {
                int a,b;
                fread(&a,sizeof(int),1,fp);
                fread(&b,sizeof(int),1,fp);
                core_t_[u][k].emplace_back(make_pair(a,b));
            }
            if (t_ < core_t_[u][k].back().second) t_ = core_t_[u][k].back().second;
        }
    }
    fclose(fp);
}

void Graph::index_baseline() {

#ifdef _LINUX_
    struct timeval t_start,t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif

    if (v_a_ == nullptr) v_a_ = new bool[n_];
    if (v_b_ == nullptr) v_b_ = new bool[n_];
    if (core_t_ == nullptr) core_t_ = new vector<vector<pair<int,int>>>[n_];
    else{
        printf("Index structure exists!\n");
        exit(1);
    }
    if (t_offset_ == nullptr) t_offset_ = new int[n_];

    core_decomposition();
    for (int u = 0; u < n_; ++u) {
        v_a_[u] = false;
        v_b_[u] = false;
        core_t_[u].resize(core_[u]+1);
    }
    memset(t_offset_,0,sizeof(int)*n_);
    compute_core_deg(0);


    for (int t_s = 0; t_s < t_; ++t_s) {

        if (t_s%100 == 0) printf("t = %d.\n",t_s);

//        for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
//            int u = edges_[i].first;
//            int v = edges_[i].second;
//
//            --nbr_cnt_[u][v];
//            if (nbr_cnt_[u][v] == 0) nbr_cnt_[u].erase(v);
//
//            --nbr_cnt_[v][u];
//            if (nbr_cnt_[v][u] == 0) nbr_cnt_[v].erase(u);
//        }
//        test_core_decomposition(t_s);

        compute_core_time_bl(t_s);
        if (t_s == t_-1) break;
        compute_core_deg(t_s);
        decremental_core_bl(t_s);
    }


#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
    printf("Running time (baseline): %lld s, %lld mins\n", t_msec/1000, t_msec/1000/60);
    if(log_f_ != nullptr) fprintf(log_f_,"Indexing time (Baseline): %lld s\n",t_msec/1000);

    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    printf("Memory usage = %ldKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024);
    if(log_f_ != nullptr) fprintf(log_f_,"Memory usage = %ldKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024);

#else
    clock_t end = clock();
    printf("Running time (baseline): %.2f s, %.2f min\n",(double)(end-start)/ CLOCKS_PER_SEC,(double)(end-start)/CLOCKS_PER_SEC/60);
#endif

    print_idx_size();
}



void Graph::compute_core_time_bl(const int &t_s) {
    vector<int> core_history;
    queue<int> q;
    int* cnt = new int[k_max_+1];
    for (int t_e = t_-1; t_e >= t_s ; --t_e) {
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e + 1]; ++i) {

            int u = edges_[i].first;
            int v = edges_[i].second;

            //        process u
            if (core_[u] <= core_[v]){
                --cd_[u][v];
                if (cd_[u][v] == 0){
                    cd_[u].erase(v);
                    if (cd_[u].size() < core_[u] && !v_a_[u]){
                        q.push(u);
                        v_a_[u] = true;
                    }
                }
            }

//            process v
            if (core_[v] <= core_[u]){
                --cd_[v][u];
                if (cd_[v][u] == 0){
                    cd_[v].erase(u);
                    if (cd_[v].size() < core_[v] && !v_a_[v]){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }
        }



        while (!q.empty()){
            int u = q.front();
            q.pop();
            v_a_[u] = false;

            int oc = core_[u];

            memset(cnt,0,sizeof(int)*(oc+1));

            for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                int t = nbr_[u][i].second;
                if (t < t_s){
                    t_offset_[u] = i+1;
                    continue;
                }
                int v = nbr_[u][i].first;
                if (t >= t_e) break;
                if (v_b_[v]) continue;
                v_b_[v] = true;
                ++cnt[core_[v] < core_[u] ? core_[v]:core_[u]];
            }

            for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b_[nbr_[u][i].first] = false;
            }

            int cd = 0;
            for (int k = oc; k >= 0 ; --k) {
                cd += cnt[k];
                if(cd >= k){
                    core_[u] = k;
                    break;
                }
            }

//            update cd_;
            cd_[u].clear();
            for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                int v = nbr_[u][i].first;
                int t = nbr_[u][i].second;
                if (t >= t_e) break;
                if (core_[v] < core_[u]) continue;
                if (cd_[u].find(v) == cd_[u].end()) cd_[u].insert(make_pair(v,1));
                else ++cd_[u][v];
            }


            //        add influenced neighbor to the queue
            for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                int v = nbr_[u][i].first;
                if (core_[u] < core_[v] && core_[v] <= oc && !v_b_[v]){
                    v_b_[v] = true;
                    cd_[v].erase(u);
                    if (!v_a_[v] && cd_[v].size() < core_[v]){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }

            for (int i = t_offset_[u]; i < nbr_[u].size(); ++i) {
                if (nbr_[u][i].second >= t_e) break;
                v_b_[nbr_[u][i].first] = false;
            }

            for (int k = oc; k > core_[u]; --k) {

                core_history.emplace_back(u);
                if (t_s == 0 || core_t_[u][k].empty() || core_t_[u][k].back().second < t_e) {
                    core_t_[u][k].emplace_back(make_pair(t_s, t_e));
                }
            }

        }

    }

    delete[] cnt;
    for (auto &i:core_history) ++core_[i];

}

void Graph::decremental_core_bl(const int &t_s) {

    queue<int> q;
    int* cnt = new int[k_max_+1];
    for (int i = edges_idx_[t_s]; i < edges_idx_[t_s + 1]; ++i) {

        int u = edges_[i].first;
        int v = edges_[i].second;

        //        process u
        if (core_[u] <= core_[v]) {
            --cd_[u][v];
            if (cd_[u][v] == 0) {
                cd_[u].erase(v);
                if (cd_[u].size() < core_[u] && !v_a_[u]) {
                    q.push(u);
                    v_a_[u] = true;
                }
            }
        }

//            process v
        if (core_[v] <= core_[u]) {
            --cd_[v][u];
            if (cd_[v][u] == 0) {
                cd_[v].erase(u);
                if (cd_[v].size() < core_[v] && !v_a_[v]) {
                    q.push(v);
                    v_a_[v] = true;
                }
            }
        }

    }
    while (!q.empty()){
        int u = q.front();
        q.pop();
        v_a_[u] = false;

        int oc = core_[u];

        memset(cnt,0,sizeof(int)*(oc+1));

        for (int i = nbr_[u].size()-1; i >= 0; --i) {
            int t = nbr_[u][i].second;
            if (t <= t_s) break;

            int v = nbr_[u][i].first;
            if (v_b_[v]) continue;
            v_b_[v] = true;
            ++cnt[core_[v] < core_[u] ? core_[v]:core_[u]];
        }

        for (int i = nbr_[u].size()-1; i >= 0; --i) {
            if (nbr_[u][i].second <= t_s) break;
            v_b_[nbr_[u][i].first] = false;
        }

        int cd = 0;
        for (int k = oc; k >= 0 ; --k) {
            cd += cnt[k];
            if(cd >= k){
                core_[u] = k;
                break;
            }
        }

//            update cd_;
        cd_[u].clear();
        for (int i = nbr_[u].size()-1; i >= 0; --i) {
            int v = nbr_[u][i].first;
            int t = nbr_[u][i].second;
            if (t <= t_s) break;
            if (core_[v] < core_[u]) continue;
            if (cd_[u].find(v) == cd_[u].end()) cd_[u].insert(make_pair(v,1));
            else ++cd_[u][v];
        }


        //        add influenced neighbor to the queue
        for (int i = nbr_[u].size()-1; i >= 0; --i) {
            if (nbr_[u][i].second <= t_s) break;
            int v = nbr_[u][i].first;
            if (core_[u] < core_[v] && core_[v] <= oc && !v_b_[v]){
                v_b_[v] = true;
                cd_[v].erase(u);
                if (!v_a_[v] && cd_[v].size() < core_[v]){
                    q.push(v);
                    v_a_[v] = true;
                }
            }
        }

        for (int i = nbr_[u].size()-1; i >= 0; --i) {
            if (nbr_[u][i].second <= t_s) break;
            v_b_[nbr_[u][i].first] = false;
        }

        for (int k = oc; k > core_[u]; --k) {
            core_t_[u][k].emplace_back(make_pair(t_s + 1, t_));
        }

    }

    delete[] cnt;

}

void Graph::online_core_decomposition(const int &t_s, const int &t_e) {

    vector<int> new_to_old;
    unordered_map<int,int> old_to_new;


    unordered_map<int,unordered_set<int>> edge_mp;
    vector<vector<int>> span_nbr;


    for (int i = edges_idx_[t_s]; i < edges_idx_[t_e+1]; ++i){
        int u = edges_[i].first;
        int v = edges_[i].second;
        if (v < u) swap(u,v);
        if (edge_mp.find(u) != edge_mp.end() && edge_mp[u].find(v) != edge_mp[u].end()) continue;

//        process u
        if (old_to_new.find(u) == old_to_new.end()){
            old_to_new.insert(make_pair(u,new_to_old.size()));
            new_to_old.emplace_back(u);
            span_nbr.emplace_back(vector<int>());
        }

//        process v
        if (old_to_new.find(v) == old_to_new.end()){
            old_to_new.insert(make_pair(v,new_to_old.size()));
            new_to_old.emplace_back(v);
            span_nbr.emplace_back(vector<int>());
        }

        span_nbr[old_to_new[u]].emplace_back(old_to_new[v]);
        span_nbr[old_to_new[v]].emplace_back(old_to_new[u]);

    }

    int span_max_deg = 0;
    for (auto &i:span_nbr){
        if (span_max_deg < i.size()) span_max_deg = i.size();
    }


    int n = new_to_old.size();
    int* span_core = new int[n];
    int* pos = new int[n];
    auto vert = new int[n];
    auto bin = new int[span_max_deg+1];
    memset(bin,0,sizeof(int)*(span_max_deg+1));

    for (int u = 0; u < n; ++u) {
        int d = span_nbr[u].size();
        span_core[u] = d;
        ++bin[d];
    }

    int offset = 0;
    for (int i = 0; i <= span_max_deg ; ++i) {
        int num = bin[i];
        bin[i] = offset;
        offset += num;
    }

    for (int u = 0; u < n; ++u) {
        pos[u] = bin[span_core[u]];
        vert[pos[u]] = u;
        bin[span_core[u]]++;
    }

    for (int i = span_max_deg; i >= 1; --i) bin[i] = bin[i - 1];
    bin[0] = 0;


    for (int i = 0; i < n; ++i) {
        int u = vert[i];

        for (auto& item : span_nbr[u]){

            if (span_core[item] > span_core[u]){
                int dv = span_core[item], pv = pos[item];
                int pw = bin[dv], w = vert[pw];
                if (item != w){
                    pos[item] = pw, vert[pv] = w;
                    pos[w] = pv, vert[pw] = item;
                }
                ++bin[dv];
                --span_core[item];
            }

        }
    }


//    span_core stores the core number in the current span

    delete[] bin;
    delete[] vert;
    delete[] pos;
    delete[] span_core;
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

void Graph::test_core_decomposition(const int &t_s) {

    if (core_ == nullptr) core_ = new int[n_];

    vector<int> old_core(n_);
    for (int u = 0; u < n_; ++u) {
        old_core[u] = core_[u];
    }

    auto pos = new int[n_];
    auto vert = new int[n_];
    auto bin = new int[max_effective_deg_+1];
    memset(bin,0,sizeof(int)*(max_effective_deg_+1));

    for (int u = 0; u < n_; ++u) {
        int d = nbr_cnt_[u].size();
        core_[u] = d;
        ++bin[d];
    }

    int offset = 0;
    for (int i = 0; i <= max_effective_deg_ ; ++i) {
        int num = bin[i];
        bin[i] = offset;
        offset += num;
    }

    for (int u = 0; u < n_; ++u) {
        pos[u] = bin[core_[u]];
        vert[pos[u]] = u;
        bin[core_[u]]++;
    }

    for (int i = max_effective_deg_; i >= 1; --i) bin[i] = bin[i - 1];
    bin[0] = 0;

    k_max_ = 0;

    for (int i = 0; i < n_; ++i) {
        int u = vert[i];

        for (auto& item : nbr_[u]){
            if (item.second < t_s) continue;
            if(v_a_[item.first]) continue;
            v_a_[item.first] = true;

            if (core_[item.first] > core_[u]){
                int dv = core_[item.first], pv = pos[item.first];
                int pw = bin[dv], w = vert[pw];
                if (item.first != w){
                    pos[item.first] = pw, vert[pv] = w;
                    pos[w] = pv, vert[pw] = item.first;
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

    delete[] bin;
    delete[] vert;
    delete[] pos;

    for (int u = 0; u < n_; ++u) {
        for (int k = old_core[u]; k > core_[u]; --k) {
            core_t_[u][k].emplace_back(make_pair(t_s,t_));
        }
    }
}

void Graph::naive_index() {

#ifdef _LINUX_
    struct timeval t_start,t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif
    long long idx_size = sizeof(int);

    for (int t_s = 0; t_s < t_; ++t_s) {
        if (t_s % 100 == 0) printf("t_s = %d\n",t_s);
        for (int t_e = t_s; t_e < t_; ++t_e) {
            online_core_decomposition(t_s,t_e);
            idx_size += n_*sizeof(int);
        }
    }

#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
    printf("Running time (Naive Index): %lld s, %lld mins\n", t_msec/1000, t_msec/1000/60);
    if(log_f_ != nullptr) fprintf(log_f_,"Indexing time (Naive Index): %lld s\n",t_msec/1000);
#else
    clock_t end = clock();
    printf("Running time (naive index): %.2f s, %.2f min\n",(double)(end-start)/ CLOCKS_PER_SEC,(double)(end-start)/CLOCKS_PER_SEC/60);
#endif

    printf("Index size: %.2f MB.\n",(float)idx_size/1024/1024);

}

void Graph::query_init() {
    if (v_a_ == nullptr) v_a_ = new bool [core_t_->size()];
    if (v_b_ == nullptr) v_b_ = new bool [core_t_->size()];
}

void Graph::query_subgraph(int u, int t_s, int t_e, int k, vector<int>& r, vector<pair<int,int>>& r_edges) {

    if (v_a_ == nullptr) query_init();

//    vector<int> r;
//    vector<pair<int,int>> r_edges;

    queue<int> q;
    if (query(u,t_s,t_e,k)){
        q.push(u);
        v_a_[u] = true;
        r.emplace_back(u);
    }
    vector<int> bm_history;

    while (!q.empty()){
        int v = q.front();
        q.pop();
        auto nbr_it = lower_bound(nbr_[v].begin(),nbr_[v].end(),make_pair(0,t_s),cmp_nbr);
        while (nbr_it != nbr_[v].end()){
            if (nbr_it->second > t_e) break;
            int w = nbr_it->first;
            ++nbr_it;
            if (v_b_[w]) continue;
            if (v_a_[w]){
                r_edges.emplace_back(make_pair(v,w));
                continue;
            }

            if (query(w,t_s,t_e,k)){
                q.push(w);
                v_a_[w] = true;
                r.emplace_back(w);
                r_edges.emplace_back(make_pair(v,w));
            }else{
                v_b_[w] = true;
                bm_history.emplace_back(w);
            }
        }

    }
    for (auto &i:bm_history) v_b_[i] = false;
    for (auto &i:r) v_a_[i] = false;

//    r is the result array



}

void Graph::init_log(const string &log_path) {
    log_f_ = fopen(log_path.c_str(),"a");
    fprintf(log_f_,"\n\n==================\n");
    time_t now = time(0);
    fprintf(log_f_,"%s\n",ctime(&now));
}

int Graph::query_all(int t_s, int t_e, int k) {
    int r = 0;
    for (int i = 0; i < n_; ++i) {
        if(query(i,t_s,t_e,k)) ++r;
    }
    return r;
}

void Graph::load(const string &path, bool timestamp_third) {
    load(path,-1,timestamp_third);
}

void Graph::naive_index_size() {

    long long idx_size = sizeof(int);

    for (int t_s = 0; t_s < t_; ++t_s) {
//        if (t_s % 100 == 0) printf("t_s = %d\n",t_s);
        for (int t_e = t_s; t_e < t_; ++t_e) {
//            online_core_decomposition(t_s,t_e);
            idx_size += n_*sizeof(int);
        }
    }

    printf("Index size (Naive Index): %.2f MB.\n",(float)idx_size/1024/1024);
    if(log_f_ != nullptr) fprintf(log_f_,"Index size (Naive Index): %.2f MB.\n",(float)idx_size/1024/1024);
}

int Graph::online_k_core(const int &t_s, const int &t_e, const int& k) {
    vector<int> new_to_old;
    unordered_map<int,int> old_to_new;


    unordered_map<int,unordered_set<int>> edge_mp;
    vector<vector<int>> span_nbr;


    for (int i = edges_idx_[t_s]; i < edges_idx_[t_e+1]; ++i){
        int u = edges_[i].first;
        int v = edges_[i].second;
        if (v < u) swap(u,v);
        if (edge_mp.find(u) != edge_mp.end() && edge_mp[u].find(v) != edge_mp[u].end()) continue;

//        process u
        if (old_to_new.find(u) == old_to_new.end()){
            old_to_new.insert(make_pair(u,new_to_old.size()));
            new_to_old.emplace_back(u);
            span_nbr.emplace_back(vector<int>());
        }

//        process v
        if (old_to_new.find(v) == old_to_new.end()){
            old_to_new.insert(make_pair(v,new_to_old.size()));
            new_to_old.emplace_back(v);
            span_nbr.emplace_back(vector<int>());
        }

        span_nbr[old_to_new[u]].emplace_back(old_to_new[v]);
        span_nbr[old_to_new[v]].emplace_back(old_to_new[u]);

    }

    int r = new_to_old.size();

    int* span_core = new int[new_to_old.size()];
    queue<int> q;
    for (int i = 0; i < span_nbr.size(); ++i) {
        span_core[i] = span_nbr[i].size();
        if (span_core[i] >= k) continue;
        q.push(i);
        new_to_old[i] = -1;
        --r;
    }

    while (!q.empty()){
        int u = q.front();
        q.pop();
        for(auto &i:span_nbr[u]){
            if (new_to_old[i] == -1) continue;
            -- span_core[i];
            if (span_core[i] >= k) continue;
            new_to_old[i] = -1;
            --r;
            q.push(i);
        }
    }
    delete[] span_core;
    return r;
}

int Graph::online_query(const int &t_s, const int &t_e, const int &k) {
    int sm;
    return online_query(t_s,t_e,k,sm);
}

int Graph::online_query(const int &t_s, const int &t_e, const int &k, int& snapshot_m) {
    vector<int> vertices;

//    reset snapshot edge number
    snapshot_m = 0;

    for (int i = edges_idx_[t_s]; i < edges_idx_[t_e+1]; ++i){
        int u = edges_[i].first;
        int v = edges_[i].second;
        if (!v_a_[u]){
            vertices.emplace_back(u);
            v_a_[u] = true;
        }
        if (!v_a_[v]){
            vertices.emplace_back(v);
            v_a_[v] = true;
        }
        query_v_[u].emplace_back(v);
        query_v_[v].emplace_back(u);
    }
    int r = vertices.size();

    queue<int> q;
    for(auto &i : vertices){
        query_deg_[i] = 0;
        vector<int> vb_history;
        for (auto &nbr: query_v_[i]){
            if (v_b_[nbr]) continue;
            v_b_[nbr] = true;
            vb_history.emplace_back(nbr);
            ++ query_deg_[i];


        }
        for(auto &nbr:vb_history) v_b_[nbr] = false;

        snapshot_m += query_deg_[i];

        if (query_deg_[i] < k){
            q.push(i);
            v_a_[i] = false;
            --r;
        }

    }

    snapshot_m /= 2;

    while (!q.empty()){
        int u = q.front();
        q.pop();
        vector<int> vb_history;
        for(auto &i:query_v_[u]){
            if (!v_a_[i] || v_b_[i]) continue;
            v_b_[i] = true;
            vb_history.emplace_back(i);
            -- query_deg_[i];
            if (query_deg_[i] < k){
                v_a_[i] = false;
                --r;
                q.push(i);
            }
        }
        for(auto &nbr:vb_history) v_b_[nbr] = false;
    }

    for (auto &u:vertices){
        query_v_[u].clear();
        v_a_[u] = false;
    }

    return r;
}

int Graph::online_span_core(const int &t_s, const int &t_e, const int &k) {
    vector<vector<pair<int,int>>> interval_edges;
    vector<pair<int,int>> intersection;


    for (int t = t_s; t <= t_e; ++t) {
        interval_edges.emplace_back(vector<pair<int,int>>());
        for (int i = edges_idx_[t]; i < edges_idx_[t+1]; ++i){
            int u = edges_[i].first;
            int v = edges_[i].second;
            if (v < u) swap(u,v);
            interval_edges.back().emplace_back(make_pair(u,v));
        }

        if (interval_edges.back().size() < k) return 0;
    }
    if (t_s == t_e){
        intersection = interval_edges.back();
    }else{
        edge_intersection(interval_edges,intersection);
        if (intersection.size() < k) return 0;
    }





    vector<int> new_to_old;
    unordered_map<int,int> old_to_new;


    unordered_map<int,unordered_set<int>> edge_mp;
    vector<vector<int>> span_nbr;

    for (auto &edge:intersection){
        int u = edge.first;
        int v = edge.second;
//        process u
        if (old_to_new.find(u) == old_to_new.end()){
            old_to_new.insert(make_pair(u,new_to_old.size()));
            new_to_old.emplace_back(u);
            span_nbr.emplace_back(vector<int>());
        }

//        process v
        if (old_to_new.find(v) == old_to_new.end()){
            old_to_new.insert(make_pair(v,new_to_old.size()));
            new_to_old.emplace_back(v);
            span_nbr.emplace_back(vector<int>());
        }

        span_nbr[old_to_new[u]].emplace_back(old_to_new[v]);
        span_nbr[old_to_new[v]].emplace_back(old_to_new[u]);

    }

    int r = new_to_old.size();

    int* span_core = new int[new_to_old.size()];
    queue<int> q;
    for (int i = 0; i < span_nbr.size(); ++i) {
        span_core[i] = span_nbr[i].size();
        if (span_core[i] >= k) continue;
        q.push(i);
        new_to_old[i] = -1;
        --r;
    }

    while (!q.empty()){
        int u = q.front();
        q.pop();
        for(auto &i:span_nbr[u]){
            if (new_to_old[i] == -1) continue;
            -- span_core[i];
            if (span_core[i] >= k) continue;
            new_to_old[i] = -1;
            --r;
            q.push(i);
        }
    }
    delete[] span_core;
    return r;
}

void Graph::edge_intersection(vector<vector<pair<int, int>>> &edges, vector<pair<int, int>> &result) {
    vector<int> pos(edges.size(),0);
    pair<int,int> cur = make_pair(-1,-1);


    for (int i = 0; i < edges.size(); ++i){
        if (pos[i] == edges[i].size()) break;
        if (cur.first == -1){
            cur.first = edges[i][pos[i]].first;
            cur.second = edges[i][pos[i]].second;
            continue;
        }
        if (edges[i][pos[i]] < cur) {
            pos[i]++;
            i--;
            continue;
        }

        if (edges[i][pos[i]] > cur) {
            cur.first = edges[i][pos[i]].first;
            cur.second = edges[i][pos[i]].second;
            i = -1;
            continue;
        }

        if (i == edges.size()-1){
            result.emplace_back(cur);
            cur.first = -1;
            pos[0] ++;
            i = -1;
        }

    }


}

int Graph::index_span_core(const int &t_s, const int &t_e, const int& k) {
    vector<vector<pair<int,int>>> interval_edges;
    vector<pair<int,int>> intersection;

    vector<int> bm_history;

    for (int t = t_s; t <= t_e; ++t) {
        interval_edges.emplace_back(vector<pair<int,int>>());
        for (int i = edges_idx_[t]; i < edges_idx_[t+1]; ++i){
            int u = edges_[i].first;
            int v = edges_[i].second;
            if (v < u) swap(u,v);

            bool edge_available = true;
            if(v_a_[u]){
                edge_available = false;
            }else if (!query(u,t,t,k)){
                v_a_[u] = true;
                bm_history.emplace_back(u);
                edge_available = false;
            }

            if(v_a_[v]){
                edge_available = false;
            }else if (!query(v,t,t,k)){
                v_a_[v] = true;
                bm_history.emplace_back(v);
                edge_available = false;
            }

            if(edge_available) interval_edges.back().emplace_back(make_pair(u,v));
        }

        if (interval_edges.back().size() < k){
            for(auto &item:bm_history) v_a_[item] = false;
            return 0;
        }
    }
    for(auto &item:bm_history) v_a_[item] = false;

    if (t_s == t_e){
        intersection = interval_edges.back();
    }else{
        edge_intersection(interval_edges,intersection);
        if (intersection.size() < k) return 0;
    }



    vector<int> new_to_old;
    unordered_map<int,int> old_to_new;


    unordered_map<int,unordered_set<int>> edge_mp;
    vector<vector<int>> span_nbr;

    for (auto &edge:intersection){
        int u = edge.first;
        int v = edge.second;
//        process u
        if (old_to_new.find(u) == old_to_new.end()){
            old_to_new.insert(make_pair(u,new_to_old.size()));
            new_to_old.emplace_back(u);
            span_nbr.emplace_back(vector<int>());
        }

//        process v
        if (old_to_new.find(v) == old_to_new.end()){
            old_to_new.insert(make_pair(v,new_to_old.size()));
            new_to_old.emplace_back(v);
            span_nbr.emplace_back(vector<int>());
        }

        span_nbr[old_to_new[u]].emplace_back(old_to_new[v]);
        span_nbr[old_to_new[v]].emplace_back(old_to_new[u]);

    }

    int r = new_to_old.size();

    int* span_core = new int[new_to_old.size()];
    queue<int> q;
    for (int i = 0; i < span_nbr.size(); ++i) {
        span_core[i] = span_nbr[i].size();
        if (span_core[i] >= k) continue;
        q.push(i);
        new_to_old[i] = -1;
        --r;
    }

    while (!q.empty()){
        int u = q.front();
        q.pop();
        for(auto &i:span_nbr[u]){
            if (new_to_old[i] == -1) continue;
            -- span_core[i];
            if (span_core[i] >= k) continue;
            new_to_old[i] = -1;
            --r;
            q.push(i);
        }
    }
    delete[] span_core;
    return r;
}


bool cmp(const pair<int,int> &a, const pair<int,int> &b){
    return a.first < b.first;
}

bool cmp_nbr(const pair<int,int> &a, const pair<int,int> &b){
    return a.second < b.second;
}

void Graph::print_index(int k) {
    cout << "index = " << endl;
    for (auto u = 0; u < n_; u++) {
        cout << "u = " << u << ": [";        
        for (auto v : core_t_[u][k]) {
            cout << "(" << v.first << "," << v.second << "),";
        }
        cout << "]" << endl;
    }
}

void Graph::time_range_kcore(long _ts, long _te, int _k) {
    cout << "Query: ts = " << _ts << ", te = " << _te << ", k = " << _k << endl;

#ifdef _LINUX_
    struct timeval t_start,t_end;
    gettimeofday(&t_start, NULL);
#else
    clock_t start = clock();
#endif

    long ts = t_old_to_new_[_ts];
    long te = t_old_to_new_[_te];
    cout << "ts = " << ts << ", te = " << te << endl;
    truncate(ts,te);

    core_t_ = new vector<vector<pair<int,int>>>[n_];    // [u1: [k1:[pair1, pair2, ...], k2:[pair1, pair2, ...], ...], u2: [], ...]
    t_offset_ = new int[n_];
    if (v_a_ == nullptr) v_a_ = new bool[n_];
    if (v_b_ == nullptr) v_b_ = new bool[n_];


    printf("starting core decomposition...\n");
//    compute core number for all edges
    core_decomposition_new();
    // cout << "core = [";
    // for (int u = 0; u < n_; u++) {
    //     cout << core_[u] << ",";
    // }
    // cout << "]" << endl;

    for (int u = 0; u < n_; ++u) {
        v_a_[u] = false;
        v_b_[u] = false;
        core_t_[u].resize(core_[u]+1);      // for each u, set size of core_t_[u] to core(u)+1
    }
    printf("k_max = %d\n",k_max_);

    printf("initialize core time.\n");
    compute_core_deg_new(ts);
    init_core_time_new(ts, te, _k);


//    init ct_deg_ and ct_cnt
    queue<int> q;
    for (int k = _k; k < _k+1; ++k) {
        for (int i = 0; i < n_; i++) {
            v_a_[i] = false;
        }
        printf("Iteration k = %d.\n",k);
        init_ct_cnt_new(k);
        for (int t_s = ts; t_s <= te; ++t_s) {
            vector<int> cand;
            for (int i = edges_idx_[t_s-1]; i < edges_idx_[t_s]; ++i) {
                int u = edges_[i].first;
                int v = edges_[i].second;

                if (invalid(u,k) || invalid(v,k)) continue;

//                process u
                if (!v_a_[u]){
                    del_nbr(u,v);
                    if (ct_cnt_[u].size()<k){
                        q.push(u);
                        v_a_[u] = true;
                    }
                }


//                process v
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

//                  compute new core time of u
                for (int i = t_offset_[u]; i < nbr_trunc_[u].size(); ++i) {

                    int t = nbr_trunc_[u][i].second;
                    if (nbr_t.size() >= k && t > ct) break;
                    if (t < t_s){
                        t_offset_[u] = i+1;
                        continue;
                    }

                    int v = nbr_trunc_[u][i].first;
//                    v is not valid or has already been visited
                    if (invalid(v,k) || v_b_[v]) continue;

//                    mark v has been visited
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
//                    sort(nbr_t.begin(),nbr_t.end());
                    new_t = nbr_t[k-1];
                }

                int old_t = core_t_[u][k].back().second;
                if (core_t_[u][k].back().first == t_s){
                    core_t_[u][k].back().second = new_t;
                }else{
                    core_t_[u][k].emplace_back(make_pair(t_s,new_t));
                }

//                compute ct_cnt_[u] and add neighbor to queue if necessary
                for (int i = t_offset_[u]; i < nbr_trunc_[u].size(); ++i) {
//                    compute ct_cnt_[u]
                    int t = nbr_trunc_[u][i].second;
                    int v = nbr_trunc_[u][i].first;
                    if (t > new_t) break;
                    if (invalid(v,k) || core_t_[v][k].back().second > new_t) continue;

                    if (new_t != t_){
                        if (ct_cnt_[u].find(v)==ct_cnt_[u].end() ){
                            ct_cnt_[u].insert(make_pair(v,1));
                        }else{
                            ++ct_cnt_[u][v];
                        }
                    }

//                    add neighbor to queue if necessary
                    if (v_a_[v]) continue;
                    if (core_t_[v][k].back().second < old_t || new_t <= core_t_[v][k].back().second) continue;
//                    del_nbr(v,u);
                    ct_cnt_[v].erase(u);
                    if (ct_cnt_[v].size() < k){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }
        }
    }





#ifdef _LINUX_
    gettimeofday(&t_end, NULL);
    long long t_msec = (t_end.tv_sec - t_start.tv_sec)*1000 + (t_end.tv_usec - t_start.tv_usec)/1000;
    printf("Running time: %lld ms, %lld s, %lld mins\n", t_msec, t_msec/1000, t_msec/1000/60);
    if(log_f_ != nullptr) fprintf(log_f_,"Indexing time: %lld s\n",t_msec/1000);


    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    printf("Memory usage = %ldKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024);
    if(log_f_ != nullptr) fprintf(log_f_,"Memory usage = %ldKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024);
#else
    clock_t end = clock();
    printf("Running time: %.2f s, %.2f min\n",(double)(end-start)/ CLOCKS_PER_SEC,(double)(end-start)/CLOCKS_PER_SEC/60);

#endif
    if(log_f_ != nullptr) fprintf(log_f_,"kmax = %d\n",k_max_);
    print_idx_size();

    //print_index(_k);
}

void Graph::truncate(long ts, long te) {
    nbr_trunc_.resize(n_);
    for (long t = ts; t <= te; t++) {
        for (int i = edges_idx_[t]; i < edges_idx_[t+1]; i++) {
            int u = edges_[i].first;
            int v = edges_[i].second;
            
            nbr_trunc_[u].emplace_back(make_pair(v,t));
            if(nbr_trunc_[u].size() > max_deg_) max_deg_trunc_ = nbr_trunc_[u].size();

            nbr_trunc_[v].emplace_back(make_pair(u,t));
            if(nbr_trunc_[v].size() > max_deg_) max_deg_trunc_ = nbr_trunc_[v].size();
        }
    }

    // cout << "nbr_trunc_ = " << endl;
    // for (int i = 0; i < nbr_trunc_.size(); i++) {
    //     cout << "u = " << i << ": [";
    //     for (int j = 0; j < nbr_trunc_[i].size(); j++) {
    //         cout << "(" << nbr_trunc_[i][j].first << "," << nbr_trunc_[i][j].second << "),";
    //     }
    //     cout << "]" << endl;
    // }

    init_nbr_cnt_trunc();

    // cout << "nbr_cnt_trunc_ = " << endl; 
    // for (int i = 0; i < n_; i++) {
    //     cout << "u = " << i << ": [";
    //     for (auto j = nbr_cnt_trunc_[i].begin(); j != nbr_cnt_trunc_[i].end(); j++) {
    //         cout << "(" << j->first << "," << j->second << "),";
    //     }
    //     cout << "]" << endl;
    // }
}

void Graph::init_nbr_cnt_trunc() {

    effective_m_ = 0;

    nbr_cnt_trunc_ = new unordered_map<int,int>[n_];


    max_effective_deg_ = 0;

    for (int u = 0; u < n_; ++u) {

        for(auto &i : nbr_trunc_[u]){
            if (nbr_cnt_trunc_[u].find(i.first) != nbr_cnt_trunc_[u].end()){
                ++nbr_cnt_trunc_[u][i.first];
            } else{
                nbr_cnt_trunc_[u].insert(make_pair(i.first,1));
            }
        }

        if(nbr_cnt_trunc_[u].size() > max_effective_deg_) max_effective_deg_ = nbr_cnt_trunc_[u].size();
        effective_m_ += nbr_cnt_trunc_[u].size();
    }

    effective_m_ /= 2;
}

//core decomposition for all edges
void Graph::core_decomposition_new() {
    if (core_ == nullptr) core_ = new int[n_];

    auto vert = new int[n_];
    auto bin = new int[max_effective_deg_+1];
    memset(bin,0,sizeof(int)*(max_effective_deg_+1));

    for (int u = 0; u < n_; ++u) {
        int d = nbr_cnt_trunc_[u].size();
        core_[u] = d;
        ++bin[d];
    }

    int offset = 0;
    for (int i = 0; i <= max_effective_deg_ ; ++i) {
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

    for (int i = 0; i < n_; ++i) {
        int u = vert[i];

        for (auto& item : nbr_trunc_[u]){
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

        for (auto& item : nbr_trunc_[u]){
            v_a_[item.first] = false;
        }

        if (core_[u] > k_max_) k_max_ = core_[u];
    }

    delete[] bin;
    delete[] vert;
    cout << "Initial decomp complete.." << endl;
}

void Graph::compute_core_deg_new(const int &t_s) {
    if (cd_ == nullptr){
        cd_ = nbr_cnt_trunc_;
    }

    for (int u = 0; u < n_; ++u) {
        cd_[u].clear();
        for (int i = nbr_trunc_[u].size()-1;i>=0;--i){
            int v = nbr_trunc_[u][i].first;
            int t = nbr_trunc_[u][i].second;
            if (t < t_s) break;

            if (core_[v] < core_[u]) continue;

            if (cd_[u].find(v) == cd_[u].end()){
                cd_[u].insert(make_pair(v,1));
            }else{
                ++cd_[u][v];
            }
        }
    }

    // cout << "cd_ = " << endl; 
    // for (int i = 0; i < n_; i++) {
    //     cout << "u = " << i << ": [";
    //     for (auto j = cd_[i].begin(); j != cd_[i].end(); j++) {
    //         cout << "(" << j->first << "," << j->second << "),";
    //     }
    //     cout << "]" << endl;
    // }
}

void Graph::init_core_time_new(long _ts, long _te, int _k) {
//    vector<int> core_history;
    int* core = new int[n_];
    for (int u = 0; u < n_; ++u) {
        core[u] = core_[u];
    }

    queue<int> q;
    int* cnt = new int[k_max_+1];
    for (int t_e = _te; t_e >= _ts; --t_e) {
        for (int i = edges_idx_[t_e]; i < edges_idx_[t_e+1]; ++i) {

            int u = edges_[i].first;
            int v = edges_[i].second;

    //        process u (delete edge)
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

//            process v (delete edge)
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
            for (int i = 0; i < nbr_trunc_[u].size(); ++i) {
                int v = nbr_trunc_[u][i].first;
                int t = nbr_trunc_[u][i].second;
                if (t >= t_e) break;
                if (v_b_[v]) continue;
                v_b_[v] = true;

                ++cnt[core[v] < core[u] ? core[v]:core[u]];
            }

            for (int i = 0; i < nbr_trunc_[u].size(); ++i) {
                if (nbr_trunc_[u][i].second >= t_e) break;
                v_b_[nbr_trunc_[u][i].first] = false;
            }

            int cd = 0;
            for (int k = oc; k >= 0 ; --k) {
                cd += cnt[k];
                if(cd >= k){
                    core[u] = k;
                    break;
                }
            }

//            update cd_;
            cd_[u].clear();
            for (int i = 0; i < nbr_trunc_[u].size(); ++i) {
                int v = nbr_trunc_[u][i].first;
                int t = nbr_trunc_[u][i].second;
                if (t >= t_e) break;
                if (core[v] < core[u]) continue;
                if (cd_[u].find(v) == cd_[u].end()) cd_[u].insert(make_pair(v,1));
                else ++cd_[u][v];
            }


    //        add influenced neighbor to the queue
            for (int i = 0; i < nbr_trunc_[u].size(); ++i) {
                if (nbr_trunc_[u][i].second >= t_e) break;
                int v = nbr_trunc_[u][i].first;
                if (core[u] < core[v] && core[v] <= oc && !v_b_[v]){
                    v_b_[v] = true;
                    cd_[v].erase(u);
                    if (!v_a_[v] && cd_[v].size() < _k){
                        q.push(v);
                        v_a_[v] = true;
                    }
                }
            }

            for (int i = 0; i < nbr_trunc_[u].size(); ++i) {
                if (nbr_trunc_[u][i].second >= t_e) break;
                v_b_[nbr_trunc_[u][i].first] = false;
            }

            core_t_[u][_k].emplace_back(make_pair(0, t_e));  // if u popped from q, then its core changed at t_e, therefore put in Index
        }
    }
    delete[] cnt;
    delete[] core;
}

void Graph::init_ct_cnt_new(int k) {

//    reuse nbr_cnt_ and avoid applying space
    ct_cnt_ = nbr_cnt_trunc_;

    for (int u = 0; u < n_; ++u) {
        if (invalid(u,k)) continue;
        t_offset_[u] = 0;
        ct_cnt_[u].clear();

        int t = core_t_[u][k].front().second;
        for (auto &i : nbr_trunc_[u]){
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

    // cout << "ct_cnt_ = " << endl; 
    // for (int i = 0; i < n_; i++) {
    //     cout << "u = " << i << ": [";
    //     for (auto j = ct_cnt_[i].begin(); j != ct_cnt_[i].end(); j++) {
    //         cout << "(" << j->first << "," << j->second << "),";
    //     }
    //     cout << "]" << endl;
    // }
}