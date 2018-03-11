/*****************************************
 *
 * Triangle Listing Algorithms
 * Author: Lijun Chang
 * Date: 05/10/2015
 * Email: ljchang@outlook.com
 *
 ******************************************
 *
 * Three categories of algorithms
 * 1) no-extra space: EdgeIterator (computing common neighbors by sort-merge join)
 * 		i) degree_oriented (optimal complexity)
 * 		ii) full_graph (baseline algorithm)
 * 		iii) least_degree_first (optimal complexity)
 *
 * 2) hash structure (O(m) space): NodeIterator (Implement min {d(u), d(v)})
 * 		i) degree_oriented (optimal complexity)
 * 		ii) full_graph (optimal complexity)
 * 		iii) least_degree_first (optimal complexity)
 *
 * 3) arrayset (O(n) space): EdgeIterator
 * 		i) degree_oriented (optimal complexity)
 * 		ii) full_graph (optimal complexity)
 *		iii) id_oriented (ad hoc algorithm)
 * 		iv) least_degree_first (optimal complexity)
 *
 * 4) compressed input (byte or nibble): EdgeIterator
 * 		i) degree_oriented (optimal complexity)
 * 		ii) full_graph (optimal complexity)
 * 		iii) id_oriented (ad hoc algorithm)
 *
 *****************************************/

#include "Timer.h"
#include "Graph.h"

using namespace std;

long long triangle_cnt;
std::vector<long long> total_time;
std::vector<long long> time1;
std::vector<long long> time2;
std::vector<long long> time3;

void print_usage() {
	printf("Usage: [1]exe [2]algorithm [3]graph-dir\n");
	printf("Algorithms:\n");
	printf("\t 1. Naive-SM\n");
	printf("\t 2. Naive-HT\n");
	printf("\t 3. TriE-BS\n");
	printf("\t 4. TriE-HT\n");
	printf("\t 5. TriC\n");
}

int main(int argc, char *argv[]) {
	if(argc < 3) {
		print_usage();
		return 0;
	}

#ifndef NDEBUG
	printf("**** Triangle Counting (Debug): %s %s *** ", argv[1], argv[2]);
#else
	printf("**** Triangle Counting (Release): %s %s *** ", argv[1], argv[2]);
#endif
	if(strcmp(argv[1], "TriE-HT") == 0) {
#ifdef _DENSE_
		printf("dense");
#else
		printf("sparse");
#endif
	}
	printf("\n");

	ui repeat_time = 5;

	for(ui i = 0;i < repeat_time;i ++) {
		Graph *graph = new Graph(argv[2]);
		graph->read_graph();
		if(strcmp(argv[1], "Naive-HT") == 0) graph->TriC_HT();
		else if(strcmp(argv[1], "TriC") == 0) graph->TriC(); //degree decreasing order
		else if(strcmp(argv[1], "TriC2") == 0) graph->TriC2(); //degree increasing order
		else if(strcmp(argv[1], "Naive-SM") == 0) graph->TriC_SM();
		else if(strcmp(argv[1], "wedge") == 0) graph->check_wedge();
		else if(strcmp(argv[1], "TriE-BS") == 0) graph->TriE_BS();
		else if(strcmp(argv[1], "TriE-HT") == 0) graph->TriE_HT();
		else print_usage();

		delete graph;
		graph = NULL;
	}

	if(total_time.size() != repeat_time) printf("WA repeat_time\n");
	if(repeat_time >= 3) {
		long long idx = 0;
		for(ui i = 1;i < repeat_time;i ++) if(total_time[i] < total_time[idx]) idx = i;
		swap(total_time[idx], total_time[repeat_time-1]);
		total_time.pop_back();
		if(!time1.empty()) {
			swap(time1[idx], time1[repeat_time-1]);
			time1.pop_back();
		}
		if(!time2.empty()) {
			swap(time2[idx], time2[repeat_time-1]);
			time2.pop_back();
		}
		if(!time3.empty()) {
			swap(time3[idx], time3[repeat_time-1]);
			time3.pop_back();
		}
		-- repeat_time;

		idx = 0;
		swap(total_time[idx], total_time[repeat_time-1]);
		total_time.pop_back();
		if(!time1.empty()) {
			swap(time1[idx], time1[repeat_time-1]);
			time1.pop_back();
		}
		if(!time2.empty()) {
			swap(time2[idx], time2[repeat_time-1]);
			time2.pop_back();
		}
		if(!time3.empty()) {
			swap(time3[idx], time3[repeat_time-1]);
			time3.pop_back();
		}
		-- repeat_time;
	}

	printf("#Triangles: %lld\n", triangle_cnt);
	long long tim = 0;
	for(ui i = 0;i < total_time.size();i ++) tim += total_time[i];
	printf("Total time: %.3lf (ms)\n", double(tim)/(repeat_time*1000));
	if(!time1.empty()) {
		tim = 0;
		for(ui i = 0;i < time1.size();i ++) tim += time1[i];
		printf("Time1: %.3lf (ms)\n", double(tim)/(repeat_time*1000));
	}
	if(!time2.empty()) {
		tim = 0;
		for(ui i = 0;i < time2.size();i ++) tim += time2[i];
		printf("Time2: %.3lf (ms)\n", double(tim)/(repeat_time*1000));
	}
	if(!time3.empty()) {
		tim = 0;
		for(ui i = 0;i < time3.size();i ++) tim += time3[i];
		printf("Time3: %.3lf (ms)\n", double(tim)/(repeat_time*1000));
	}

	// printf("\t*** Finished triangulation\n");
	printf("\n");

	return 0;
}
