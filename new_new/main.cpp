#include "newGraph.h"

#include <iostream>
#include <random>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {
    string graph_path(argv[1]);

    long ts = stol(argv[2]);
    long te = stol(argv[3]);
    int k = stoi(argv[4]);
    string version = argv[5];
    int write_res = stoi(argv[6]);
    int threads = stoi(argv[7]);

    auto *g = new Graph();
    g->load(graph_path, ts, te, k, threads, write_res);

    if (version == "v1") {
        g->baseline_serial();
    }
    else if (version == "v2") {
        g->decrement_parallel();
    }
    else if (version == "v3") {
        g->queue_parallel();
    }
    else if (version == "v4") {
        g->vcentric_parallel();
    }

    delete g;
    return 0;
}