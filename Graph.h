#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "HashSet.h"

class Graph {
private:
	std::string dir; //input graph directory
	ui n, m; //number of nodes and edges of the graph

	ui *pstart; //offset of neighbors of nodes
	ui *edges; //adjacent ids of edges

#ifdef _HASHSET_
	HashSet **hs;
	//unordered_set<int> **hs;
#else
	hash_set<ui> **hs;
	hash_map<ui, ui> **ht;
#endif

public:
	Graph(const char *_dir) ;
	~Graph() ;

	void read_graph() ;

	void check_wedge() ;

	void ni_oriented() ;
	void TriC_HT() ;
	void TriE_HT();
	void TriE_BS();
	void TriC() ;
	void TriC2() ;
	void TriC_SM() ;
	void ei_full_co() ;
	void ei_oriented_v() ;
	void ei_roriented_v() ;
	void ei_oriented_c() ;

	void destroy_hash() ;

private:
	inline ui binary_search(const ui *array, ui e, ui val) ;
	inline void cross_link(ui *reverse) ;
	void build_degree_oriented_graph() ;
	void build_degree_oriented_graph_reverse() ;
	void build_least_oriented_graph() ;

	void increment(ui u, ui v, ui *tri_cnt, ui *deg, ui *pend) ;

	void build_hashset() ;
	void build_hashtable(ui *pend) ;
	int find_hashset(ui a, ui b) ;
	ui find_hashtable(ui u, ui v) ;
	void increment_hashtable(ui u, ui v) ;
};

#endif
