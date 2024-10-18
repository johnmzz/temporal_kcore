#include "newGraph.h"
#include <iostream>
#include <random>
#include <cmath>
#include <sys/resource.h> // Required for getrusage

#ifdef __linux__
#include <unistd.h>
#endif

#ifdef __APPLE__
#include <mach/mach.h>
#endif

using namespace std;

// Function to get the current memory usage (RSS)
long getCurrentMemoryUsage() {
#ifdef __linux__
    long rss = 0L;
    std::ifstream statm("/proc/self/statm");
    if (statm.is_open()) {
        long pages;
        statm >> pages; // Read total program size (pages)
        statm >> pages; // Read resident set size (pages)
        rss = pages * sysconf(_SC_PAGESIZE) / 1024L; // Convert pages to kilobytes
    }
    return rss;
#elif defined(__APPLE__)
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) != KERN_SUCCESS) {
        return 0L; // Failed to get memory usage
    }
    return info.resident_size / 1024L; // Convert bytes to kilobytes
#else
    return 0L; // Unsupported platform
#endif
}

// Function to get the maximum resident set size (Max RSS)
long getMaxMemoryUsage() {
#ifdef __linux__
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // Max RSS in kilobytes
#elif defined(__APPLE__)
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) != KERN_SUCCESS) {
        return 0L; // Failed to get memory usage
    }
    return info.resident_size / 1024L; // Convert bytes to kilobytes
#else
    return 0L; // Unsupported platform
#endif
}

void printMemoryUsage() {
    long currentMemoryUsage = getCurrentMemoryUsage();
    long maxMemoryUsage = getMaxMemoryUsage();
    
    std::cout << "Current Memory Usage (RSS): " << currentMemoryUsage << " KB, "
              << currentMemoryUsage / 1024.0 << " MB, "
              << currentMemoryUsage / (1024.0 * 1024.0) << " GB" << std::endl;
    
    std::cout << "Max Memory Usage (Max RSS): " << maxMemoryUsage << " KB, "
              << maxMemoryUsage / 1024.0 << " MB, "
              << maxMemoryUsage / (1024.0 * 1024.0) << " GB" << std::endl;
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <graph_path> <ts> <te> <k> <algorithm>" << std::endl;
        return 1;
    }

    string graph_path(argv[1]);
    long ts = stol(argv[2]);
    long te = stol(argv[3]);
    int k = stoi(argv[4]);
    string version = argv[5];
    int write_res = stoi(argv[6]);
    string algorithm = argv[7];

    // Step 1: Load the graph
    auto *g = new Graph();
    g->load(graph_path, ts, te, k, version, write_res);

    // Perform operations and measure memory usage
    if (version == "v") {
        g->core_time();
    } else if (version == "e") {
        g->e_core_time();
    }

    if (algorithm == "baseline") {
        g->baseline();
    } else if (algorithm == "advanced") {
        g->compute_sets();
    }
    printMemoryUsage();

    // Delete the graph
    delete g;
    return 0;
}
