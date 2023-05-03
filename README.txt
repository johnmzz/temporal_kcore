READ ME.txt
Contains code and scripts relating to submission "On Querying Hisotrical K-Cores"

------------------
IMPORTANT: For different temporal graph dataset formats. 
------------------
by default, the expected datafile contains three columns <u> <v> <timestamp>
e.g. 1 2 1082040961
use -idx if data follows above format with weight (i.e. Youtube, DBLP) 
$ ./run -idx ../dataset/CollegeMsg.txt ../idx/CollegeMsg_idx ./results-rest/new-log-CollegeMsg-idx.txt_idx

for other datasets, an additional column for the edge weight is common <u> <v> <weight> <timestamp>
e.g. 1	2	1 410227261  
use -idx-4col if data follows above format with weight (i.e. Youtube, DBLP) 
$ ./run -idx-4col ../dataset/youtube-u-growth/out.youtube-u-growth ../idx/youtube_idx_d ./results-rest/nnew-log-youtube-idx.txt

------------------
Index Contruction:
Runs HC*-Construct* OR HC*-Construct
------------------
Arguments: 
./run -idx <path_to_graph> <path_to_index> <log_file> (<number_edges>) OR
./run -idx-bl <path_to_graph> <path_to_index> <log_file> (<number_edges>)
0) ./run
1) flags 	
	if '-idx', optimized HC*-Construct* algorithm 
	if '-idx-bl' HC*-Construct algorithm
2) path to the graph
3) path of index destination 
4) log filename (appends to file if exists, else new file is created)
5) optional arguement (#edges) defaults to all edges in graph

e.g.
$ ./run -idx-bl ../dataset/dblp_coauthor/out.dblp_coauthor ../idx/dblp_idx_scale dblp_idx_log 29487744

------------------
(If index file exists) Querying Index: 
Runs Online-Query AND HC*-Query
------------------
generates a specified number of queries (see arg 6)
user controls the size of the time window (see arg 5)
the code will go through all preset k-size (20%, 40%, 60%, 80%)
for each t and k, log file reports the average time from all queries (both Online-Query and HC*-Query)

Arguments:
./run -q <path_to_graph> <path_to_index> <log_file> <number_edges>
0) ./run
1) flag '-q'
2) path to the graph
3) path of index destination 
4) log filename (appends to file if exists, else creates new file)
5) integer value (1-9) that controls size of the window (e.g. 1 = t_max*0.1)
6) number of queries to generate/run

e.g.
$ ./run -q ../dataset/CollegeMsg.txt ../idx/CollegeMsg_idx college_log 6 1000

	
