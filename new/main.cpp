#include "myGraph.h"

#include <iostream>
#include <random>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {
    // g++ -O3 main.cpp myGraph.cpp -o run -std=c++11
    // clang++ -O3 -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp -std=c++11 main.cpp myGraph.cpp -o run
    // g++ -O3 -fopenmp -std=c++11 main.cpp myGraph.cpp -o run
    // ./run ../data/graph.txt idx-log-CollegeMsg.txt 1 8 2 1
    // ./run ../data/CollegeMsg.txt idx-log-CollegeMsg.txt 1082040961 1098777142 2 1
    // ./run ../data/CollegeMsg_day.txt idx-log-CollegeMsg.txt 0 193 2 1
    string graph_path(argv[1]);
    string log_path(argv[2]);

    long ts = stol(argv[3]);
    long te = stol(argv[4]);
    int k = stoi(argv[5]);
    int threads = stoi(argv[6]);

    auto *g = new Graph();
    g->init_log(log_path);
    g->load(graph_path);

    // g->time_range_kcore(ts, te, k);

    // g->time_range_kcore_parallel(ts, te, k, threads);

    g->time_range_kcore_v2(ts, te, k, threads);
    // g->print_graph();

    delete g;
    return 0;
}

