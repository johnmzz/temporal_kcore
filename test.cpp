#include <iostream>
#include <vector>
#include <omp.h>

using namespace std;

int main() {
    int total = 0;
    vector<vector<int>> table = vector<vector<int>>(5, vector<int>(5));
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            table[i][j] = i;
        }
    }
    for (int i = 0; i < 5; i++) {
        cout << "[";
        for (int j = 0; j < 5; j++) {
            cout << table[i][j] << ",";
        }
        cout << "]\n";
    } cout << endl;

    int* table2 = new int[10000];
    for (int i = 0; i < 10000; i++) {
        table2[i] = 100;
    }
    cout << table2[50];

    vector<int> table3 = vector<int>(10000, -1);
    cout << table3[50];


    // #pragma omp parallel for num_threads(10)
    // for (int i = 0; i < 4; i++) {
    //     for (int j = 0; j < 5; j++) {
    //         table[i][j] += table[i+1][j];
    //     }
    //     total += omp_get_thread_num();
    //     // cout << "Hello OpenMP from thread " << omp_get_thread_num() << " out of " << omp_get_num_threads() << " threads." << endl;
    // }

    for (int i = 0; i < 5; i++) {
        cout << "[";
        for (int j = 0; j < 5; j++) {
            cout << table[i][j] << ",";
        }
        cout << "]\n";
    }
    cout << "total = " << total << endl;

    struct rusage rUsage;
    getrusage(RUSAGE_SELF, &rUsage);
    long ms = rUsage.ru_maxrss;
    printf("Memory usage = %ld B, %.2fKB, %.2fMB, %.2fGB\n",ms,(float)ms/1024,(float)ms/1024/1024,(float)ms/1024/1024/1024);
    return 0;
}