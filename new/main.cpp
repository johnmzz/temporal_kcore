#include "myGraph.h"

#include <iostream>
#include <random>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {
    // g++ -O3 main.cpp myGraph.cpp -o run -std=c++11
    // clang++ -O3 -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp -std=c++11 main.cpp myGraph.cpp -o run
    // g++ -O3 -fopenmp -std=c++11 main.cpp myGraph.cpp -o run
    // ./run ../data/graph.txt 1 8 2 1
    // ./run ../data/CollegeMsg.txt 1082040961 1098777142 2 v4 1
    // ./run ../data/CollegeMsg_day.txt 0 193 2 v4 1
    // ./run ../data/email_day.txt 0 803 2 1
    // ./run ../data/mathoverflow_day.txt 0 2350 2 1
    // ./run ../data/askUbuntu_day.txt 0 2613 2 1
    // ./run ../data/superuser_day.txt 0 2773 2 1
    // ./run ../data/WikiTalk_day.txt 0 2320 2 1
    // ./run ../data/flickr.txt 1162422000 1164927600 2 v3 1
    string graph_path(argv[1]);

    long ts = stol(argv[2]);
    long te = stol(argv[3]);
    int k = stoi(argv[4]);
    string version = argv[5];
    int threads = stoi(argv[6]);

    auto *g = new Graph();
    g->load(graph_path);

    if (version == "v3") {
        g->time_range_kcore_v3(ts, te, k, threads);
    }
    else if (version == "v4") {
        g->time_range_kcore_v4(ts, te, k, threads);
    }
    else if (version == "v1") {
        g->time_range_kcore(ts, te, k);
    }

    delete g;
    return 0;
}

















