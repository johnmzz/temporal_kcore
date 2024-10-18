source code for the paper: Accelerating K-Core Computation in Temporal Graphs
# prerequisites
 - gcc 4.8.5 or higher
# build the code
``
g++ -O3 main.cpp newGraph.cpp -o run -std=c++14
``
# run the demo
``
./run [data_file] [start_time] [end_time] [k] [algoritm]
``

For algorithms, input "advanced" or "baseline"

For example:
``
./run ../data/CollegeMsg.txt 1082040961 1098777142 2 advanced
``

Data file format (txt file, each line should look like this):
``
u v timestamp
``

Please make sure all temporal graph files are sorted by increasing order of timestamps.
