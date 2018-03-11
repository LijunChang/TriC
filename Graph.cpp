#include "Graph.h"

using namespace std;

extern long long triangle_cnt;
extern std::vector<long long> total_time;
extern std::vector<long long> time1;
extern std::vector<long long> time2;
extern std::vector<long long> time3;

Graph::Graph(const char *_dir) {
	dir = string(_dir);

	n = m = 0;

	pstart = NULL;
	edges = NULL;

	hs = NULL;
#ifndef _HASHSET_
	ht = NULL;
#endif
}

Graph::~Graph() {
	if(pstart != NULL) {
		delete[] pstart;
		pstart = NULL;
	}
	if(edges != NULL) {
		delete[] edges;
		edges = NULL;
	}

	if(hs != NULL) {
		for(ui i = 0;i < n;i ++) if(hs[i] != NULL) {
			delete hs[i];
			hs[i] = NULL;
		}
		delete[] hs;
		hs = NULL;
	}

#ifndef _HASHSET_
	if(ht != NULL) {
		for(ui i = 0;i < n;i ++) if(ht[i] != NULL) {
			delete ht[i];
			ht[i] = NULL;
		}
		delete[] ht;
		ht = NULL;
	}
#endif
}

void Graph::destroy_hash() {
	if(hs != NULL) {
		for(ui i = 0;i < n;i ++) if(hs[i] != NULL) {
			delete hs[i];
			hs[i] = NULL;
		}
		delete[] hs;
		hs = NULL;
	}

#ifndef _HASHSET_
	if(ht != NULL) {
		for(ui i = 0;i < n;i ++) if(ht[i] != NULL) {
			delete ht[i];
			ht[i] = NULL;
		}
		delete[] ht;
		ht = NULL;
	}
#endif
}

void Graph::read_graph() {
	FILE *f = open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != (int)sizeof(int)) {
		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	printf("\tn = %u; m = %u\n", n, m/2);

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, f);

#ifndef NDEBUG
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	if(sum != m) printf("WA input graph\n");
#endif

	fclose(f);

	f = open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(pstart == NULL) pstart = new ui[n+1];
	if(edges == NULL) edges = new ui[m];

	srandom(0);

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {

		if(degree[i] > 0) {
			fread(edges+pstart[i], sizeof(int), degree[i], f);
			for(ui j = 0;j < degree[i];j ++) if(edges[pstart[i]+j] == i) {
				edges[pstart[i]+j] = edges[pstart[i]+degree[i]-1];
				-- degree[i];
			}

		}
		pstart[i+1] = pstart[i] + degree[i];
		random_shuffle(edges+pstart[i], edges+pstart[i+1]);
	}

	fclose(f);

	//FILE *fout = open_file("log.bin", "wb");
	//fwrite(edges, sizeof(int), m, fout);
	//fclose(fout);

	int self_loop = 0, not_sorted = 0;

	for(ui i = 0;i < n;i ++) {
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(edges[j] == i) self_loop = 1;
			if(j > pstart[i]&&edges[j] <= edges[j-1]) not_sorted = 1;
		}
	}

	if(self_loop) printf("!!! Self_loop\n");
	//if(not_sorted) printf("!!! Not_sorted\n");

	delete[] degree;
}

void Graph::build_degree_oriented_graph() {
	ui *degree = new ui[n+1];
	memset(degree, 0, sizeof(ui)*(n+1));

	for(ui i = 0;i < n;i ++) {
		ui d = pstart[i+1] - pstart[i];
		++ degree[d];
	}
	for(ui i = 1;i < n;i ++) degree[i] += degree[i-1];

	ui *oid = new ui[n];
	for(ui i = 0;i < n;i ++) {
		ui d = pstart[i+1] - pstart[i];
		oid[i] = degree[d];
		-- degree[d];
	}

	ui *new_pstart = degree;

	ui cid = 0;
	for(ui i = 0;i < n;i ++) {
		new_pstart[i] = cid;
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(oid[edges[j]] > oid[i]) {
				edges[cid] = edges[j];
				++ cid;
			}
		}
	}
	new_pstart[n] = cid;
	for(ui i = 0;i <= n;i ++) pstart[i] = new_pstart[i];

	delete[] degree;
	delete[] oid;
}

void Graph::build_degree_oriented_graph_reverse() {
	ui *degree = new ui[n+1];
	memset(degree, 0, sizeof(ui)*(n+1));

	for(ui i = 0;i < n;i ++) {
		ui d = pstart[i+1] - pstart[i];
		++ degree[d];
	}
	for(ui i = 1;i < n;i ++) degree[i] += degree[i-1];

	ui *oid = new ui[n];
	for(ui i = 0;i < n;i ++) {
		ui d = pstart[i+1] - pstart[i];
		oid[i] = degree[d];
		-- degree[d];
	}

	ui *new_pstart = degree;

	ui cid = 0;
	for(ui i = 0;i < n;i ++) {
		new_pstart[i] = cid;
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(oid[edges[j]] < oid[i]) {
				edges[cid] = edges[j];
				++ cid;
			}
		}
	}
	new_pstart[n] = cid;
	for(ui i = 0;i <= n;i ++) pstart[i] = new_pstart[i];

	delete[] degree;
	delete[] oid;
}

void Graph::build_least_oriented_graph() {
	ui *degree = new ui[n];
	ui *start = new ui[n+1];
	memset(start, 0, sizeof(ui)*(n+1));

	for(ui i = 0;i < n;i ++) {
		ui d = pstart[i+1] - pstart[i];
		degree[i] = d;
		++ start[d];
	}
	for(ui i = 1;i <= n;i ++) start[i] += start[i-1];

	ui *oid = new ui[n];
	ui *ver = new ui[n];
	for(ui i = 0;i < n;i ++) {
		ui d = pstart[i+1] - pstart[i];
		oid[i] = -- start[d];
		ver[start[d]] = i;
	}

	for(ui i = 0;i < n;i ++) {
		ui u = ver[i];
		start[degree[u]] ++;
		if(degree[u] > 0) start[degree[u]-1] = start[degree[u]];

		for(ui j = pstart[u];j < pstart[u+1];j ++) {
			ui v = edges[j];
			if(oid[v] <= i) continue;

			if(start[degree[v]] != oid[v]) {
#ifdef _DEBUG_
				ui u = ver[start[degree[v]]];
#endif
				swap(ver[start[degree[v]]], ver[oid[v]]);
				oid[ver[oid[v]]] = oid[v];
				oid[v] = start[degree[v]];
#ifdef _DEBUG_
				if(ver[oid[u]] != u||ver[oid[v]] != v) printf("WA\n");
#endif
			}

			++ start[degree[v]];
			-- degree[v];
		}
	}

	ui *new_pstart = start;

	ui cid = 0;
	for(ui i = 0;i < n;i ++) {
		new_pstart[i] = cid;
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(oid[edges[j]] > oid[i]) {
				edges[cid] = edges[j];
				++ cid;
			}
		}
	}

	new_pstart[n] = cid;
	for(ui i = 0;i <= n;i ++) pstart[i] = new_pstart[i];

	delete[] degree;
	delete[] start;
	delete[] oid;
	delete[] ver;
}

/*
void Graph::ni_oriented() {
#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#else
	int start, end1, end;
	start = clock();
#endif

	build_hashset();

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#else
	end1 = clock();
#endif

	if(ig == degree_oriented) build_degree_oriented_graph();
	else if(ig == least_oriented) build_least_oriented_graph();
	else printf("WA inputgraph\n");

	//printf("Finished orientation\n");

#ifdef _LINUX_
	struct timeval end2;
	gettimeofday(&end2, NULL);

	long long mtime2, seconds2, useconds2;
	seconds2 = end2.tv_sec - end1.tv_sec;
	useconds2 = end2.tv_usec - end1.tv_usec;
	mtime2 = seconds2*1000000 + useconds2;
#else
	int end2 = clock();
#endif

	long long res = 0;
	for(ui i = 0;i < n;i ++) {
		//if(i%1000 == 0) printf("%d\n", i);
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			for(ui k = j+1;k < pstart[i+1];k ++) {
				if(find_hashset(edges[j], edges[k])) ++ res;
			}
		}
	}

	printf("\t#Triangles by NI-Oriented: %lld\n", res);

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	//printf("\tTriangle time: %lld\n", mtime-mtime1-mtime2);
	printf("Hash time:\t%lld\nOrient time:\t%lld\nTriangle time:\t%lld\nTotal time:\t%lld\n", mtime1, mtime2, mtime-mtime1-mtime2, mtime);
#else
	end = clock();

	printf("Hash time: %d, Orientation time: %d, Triangle time: %d\n", end1-start, end2-end1, end-end2);
#endif
}
*/

void Graph::TriC_HT() {
	Timer t;

	build_hashset() ;

	long long mtime1 = t.elapsed();

	long long res = 0;
	ui *tri_cnt = new ui[m];

	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pstart[i+1];j ++) {
		ui &cnt = tri_cnt[j] = 0;
		ui v = edges[j];
		ui x = i, y = v;
		if(pstart[v+1] - pstart[v] < pstart[i+1] - pstart[i]) {
			x = v;
			y = i;
		}

		for(ui k = pstart[x];k < pstart[x+1];k ++) {
			if(find_hashset(y, edges[k])) ++ cnt;
		}

		res += cnt;
	}

	triangle_cnt = res;
	total_time.pb(t.elapsed());

	//printf("#Triangles by TriC-HT: %lld\nHash_Time: %lld, Total_Time: %lld\n", res, mtime1, t.elapsed());

	delete[] tri_cnt;
}

void Graph::check_wedge() {
	ui *deg = new ui[n];
	ui *degp = new ui[n];
	ui *degm = new ui[n];

	for(ui i = 0;i < n;i ++) deg[i] = pstart[i+1]-pstart[i];

	for(ui i = 0;i < n;i ++) {
		degp[i] = degm[i] = 0;
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			ui v = edges[j];
			if(deg[v] > deg[i]||(deg[v] == deg[i]&&v < i)) ++ degp[i];
			if(v < i) ++ degm[i];
		}
	}

	long long trism = 0, trihs = 0, trio = 0, trif = 0;
	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pstart[i+1];j ++) {
		ui v = edges[j];
		trism += deg[i] + deg[v];

		if(deg[i] < deg[v]) trihs += deg[i];
		else trihs += deg[v];

		if(deg[v] > deg[i]||(deg[v] == deg[i]&&v < i)) trio += degp[v];

		if(deg[v] < deg[i]||(deg[v] == deg[i]&&v < i)) trif += degm[v];
	}

	printf("trism: %lld\n", trism);
	printf("trihs: %lld\n", trihs);
	printf("trioc: %lld\n", trio);
	printf("trifc: %lld\n", trif);

	delete[] deg;
	delete[] degp;
	delete[] degm;
}

void Graph::TriC_SM() {
	Timer t;

	for(ui i = 0;i < n;i ++) sort(edges+pstart[i], edges+pstart[i+1]);

	long long mtime1 = t.elapsed();

	long long res = 0;
	ui *tri_cnt = new ui[m];

	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pstart[i+1];j ++) {
		ui &cnt = tri_cnt[j] = 0;
		ui v = edges[j];

		ui ii = pstart[i], jj = pstart[v];
		while(ii < pstart[i+1]&&jj < pstart[v+1]) {
			if(edges[ii] < edges[jj]) ++ ii;
			else if(edges[ii] > edges[jj]) ++ jj;
			else {
				++ cnt;
				++ ii; ++ jj;
			}
		}
		res += cnt;
	}

	triangle_cnt = res;
	total_time.pb(t.elapsed());

	//long long mtime2 = t.elapsed();
	//printf("#Triangles by TriC-SM: %lld\nTri_Time: %lld, Total_Time: %lld\n", res, mtime2-mtime1, mtime2);

	delete[] tri_cnt;
}

void Graph::ei_full_co() {
#ifdef _LINUX_
	struct timeval start;
	gettimeofday(&start, NULL);
#else
	int start = clock();
#endif

	long long res = 0;

	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pstart[i+1];j ++) {
		if(edges[j] < i) continue;

		ui v = edges[j];

		int ii = pstart[i+1]-1, jj = pstart[v+1]-1;
		while(ii >= pstart[i]&&jj >= pstart[v]&&edges[ii] > v&&edges[jj] > v) {
			if(edges[ii] < edges[jj]) -- jj;
			else if(edges[ii] > edges[jj]) -- ii;
			else {
				-- ii; -- jj;
				++ res;
			}
		}
	}

	printf("\t#Triangles by EI-FULL-C: %lld\n", res);

#ifdef _LINUX_
	struct timeval end1;
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;

	printf("Total time:\t%lld\n", mtime1);
#else
	int end1 = clock();
#endif
}

ui Graph::binary_search(const ui *array, ui e, ui val) {
	assert(e > 0);
	ui b = 0;
	-- e;
#ifndef NDEBUG
	if(array[e] < val||array[0] > val) {
		printf("Not found in binary_search!\n");
		return e+1;
	}
#endif
	while(b < e) {
		ui mid = b + (e-b)/2;
		if(array[mid] >= val) e = mid;
		else b = mid+1;
	}
	return e;
}

void Graph::cross_link(ui *reverse) {
	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pstart[i+1];j ++) {
		ui v = edges[j];
		ui d1 = pstart[i+1] - pstart[i], d2 = pstart[v+1] - pstart[v];

		if(d1 < d2||(d1 == d2&&i < v)) continue;

		ui r_id = binary_search(edges+pstart[v], pstart[v+1]-pstart[v], i) + pstart[v];
#ifndef NDEBUG
		if(r_id < pstart[v]||r_id >= pstart[v+1]||i != edges[r_id]) printf("??? WA in cross_link\n");
#endif

		reverse[j] = r_id;
		reverse[r_id] = j;
	}
}

void Graph::increment(ui u, ui v, ui *tri_cnt, ui *deg, ui *pend) {
	if(deg[v] < deg[u]||(deg[v] == deg[u]&&v < u)) {
		ui idx = pstart[u] + binary_search(edges+pstart[u], pend[u] - pstart[u], v);
		++ tri_cnt[idx];
	}
	else {
		ui idx = pend[u] + binary_search(edges+pend[u], pstart[u+1]-pend[u], v);
		++ tri_cnt[idx];
	}
}

void Graph::TriE_BS() {
	Timer t;
	ui *deg = new ui[n];
	for(ui i = 0;i < n;i ++) deg[i] = pstart[i+1]-pstart[i];

	ui *pend = new ui[n];
	for(ui i = 0;i < n;i ++) {
		ui &end = pend[i] = pstart[i];
		ui j = pstart[i+1];
		while(true) {
			while(end < j&&(deg[edges[end]] < deg[i]||(deg[edges[end]]==deg[i]&&edges[end]<i))) ++ end;
			while(j > end&&(deg[edges[j-1]] > deg[i]||(deg[edges[j-1]]==deg[i]&&edges[j-1]>i))) -- j;
			if(end >= j) break;
			swap(edges[end], edges[j-1]);
		}
		sort(edges+pstart[i], edges+pend[i]);
		sort(edges+pend[i], edges+pstart[i+1]);
	}

	long long mtime1 = t.elapsed();

	ui *tri_cnt = new ui[m];
	memset(tri_cnt, 0, sizeof(ui)*m);

	char *adj = new char[n];
	memset(adj, 0, sizeof(char)*n);

	long long res = 0;
	for(ui u = 0;u < n;u ++) {
		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 1;

		for(ui j = pstart[u];j < pend[u];j ++) {
			int v = edges[j];

			for(ui k = pstart[v];k < pend[v];k ++) {
				ui w = edges[k];
				if(adj[w]) {
					increment(u, v, tri_cnt, deg, pend);
					increment(v, u, tri_cnt, deg, pend);
					increment(v, w, tri_cnt, deg, pend);
					increment(w, v, tri_cnt, deg, pend);
					increment(u, w, tri_cnt, deg, pend);
					increment(w, u, tri_cnt, deg, pend);
					res += 6;
				}
			}
		}

		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}

	triangle_cnt = res;
	total_time.pb(t.elapsed());

	//long long mtime3 = t.elapsed();
	//printf("#Triangles by TriE-BS: %lld\nTotal_Time: %lld\n", res, mtime3);
	//printf("Sort_Time: %lld, Tri_Time: %lld\n", mtime1, mtime3-mtime1);

	delete[] tri_cnt;
	delete[] pend;
	delete[] deg;
	delete[] adj;
}

void Graph::TriE_HT() {
	Timer t;

	ui *deg = new ui[n];
	for(ui i = 0;i < n;i ++) deg[i] = pstart[i+1]-pstart[i];

	ui *pend = new ui[n];
	for(ui i = 0;i < n;i ++) {
		ui &end = pend[i] = pstart[i];
		ui j = pstart[i+1];
		while(true) {
			while(end < j&&(deg[edges[end]] < deg[i]||(deg[edges[end]]==deg[i]&&edges[end]<i))) ++ end;
			while(j > end&&(deg[edges[j-1]] > deg[i]||(deg[edges[j-1]]==deg[i]&&edges[j-1]>i))) -- j;
			if(end >= j) break;
			swap(edges[end], edges[j-1]);
		}
	}

	long long mtime1 = t.elapsed();

	build_hashtable(pend);

	long long mtime2 = t.elapsed();

	ui *tri_cnt = new ui[m];
	memset(tri_cnt, 0, sizeof(ui)*m);

	char *adj = new char[n];
	memset(adj, 0, sizeof(char)*n);

	//long long tmp = 0;
	long long res = 0;
	for(ui u = 0;u < n;u ++) {
		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 1;

		for(ui j = pstart[u];j < pend[u];j ++) {
			int v = edges[j];

			for(ui k = pstart[v];k < pend[v];k ++) {
				ui w = edges[k];
				if(adj[w]) {
					increment_hashtable(u,v);
					increment_hashtable(v,w);
					increment_hashtable(u,w);

					res += 6;
				}
			}
		}

		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}

	long long mtime4 = t.elapsed();

	for(ui i = 0;i < n;i ++) {
		for(ui j = pstart[i];j < pend[i];j ++) tri_cnt[j] = find_hashtable(i, edges[j]);
		for(ui j = pend[i];j < pstart[i+1];j ++) tri_cnt[j] = find_hashtable(edges[j], i);
	}

	long long mtime3 = t.elapsed();

	triangle_cnt = res;
	total_time.pb(mtime3);
	time1.pb(mtime2);
	time2.pb(mtime4-mtime2);
	time3.pb(mtime3-mtime4);

	//printf("#Triangles by TriE-HT: %lld\nTotal_Time: %lld\n", res, mtime3);
	//printf("Organization_Time: %lld, Hash_Time: %lld, Tri_Time: %lld, cross time: %lld\n", mtime1, mtime2-mtime1, mtime4-mtime2, mtime3-mtime4);

	delete[] tri_cnt;
	delete[] pend;
	delete[] deg;
	delete[] adj;
}

void Graph::TriC() {
	Timer t;

	ui *adj = new ui[n];
	ui *deg = adj;
	for(ui i = 0;i < n;i ++) deg[i] = pstart[i+1]-pstart[i];

	ui *pend = new ui[n];

	for(ui i = 0;i < n;i ++) {
		ui &end = pend[i] = pstart[i];
		ui j = pstart[i+1];
		while(true) {
			while(end < j&&(deg[edges[end]] < deg[i]||(deg[edges[end]]==deg[i]&&edges[end]<i))) ++ end;
			while(j > end&&(deg[edges[j-1]] > deg[i]||(deg[edges[j-1]]==deg[i]&&edges[j-1]>i))) -- j;
			if(end >= j) break;
			swap(edges[end], edges[j-1]);
		}
		sort(edges+pend[i], edges+pstart[i+1]);
	}

	long long mtime1 = t.elapsed();

	ui *tri_cnt = new ui[m];
	memset(tri_cnt, 0, sizeof(ui)*m);

	memset(adj, 0, sizeof(ui)*n);

	long long res = 0;
	for(ui u = 0;u < n;u ++) {
		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = j+1;

		for(ui j = pstart[u];j < pend[u];j ++) {
			int v = edges[j];

			for(ui k = pstart[v];k < pend[v];k ++) {
				//++ tmp;
				if(adj[edges[k]]) {
					++ tri_cnt[j];
					++ tri_cnt[k];
					++ tri_cnt[adj[edges[k]]-1];
					res += 6;
				}
			}
		}

		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}

	long long mtime2 = t.elapsed();

	for(ui u = 0;u < n;u ++) {
		for(ui j = pstart[u];j < pend[u];j ++) {
			ui v = edges[j];

			ui r_id = binary_search(edges+pend[v], pstart[v+1]-pend[v], u) + pend[v];
			tri_cnt[r_id] = tri_cnt[j];
		}
	}

	long long mtime3 = t.elapsed();
	triangle_cnt = res;
	total_time.pb(mtime3);
	time1.pb(mtime1);
	time2.pb(mtime2-mtime1);
	time3.pb(mtime3-mtime2);

	//printf("#Triangles by TriC: %lld\nTotal_Time: %lld\n", res, mtime3);
	//printf("Sort_Time: %lld, CL_Time: %lld, Tri_Time: %lld\n", mtime1, mtime3-mtime2, mtime2-mtime1);

	delete[] tri_cnt;
	delete[] pend;
	delete[] adj;
}

void Graph::TriC2() {
	Timer t;

	ui *adj = new ui[n];
	ui *deg = adj;
	for(ui i = 0;i < n;i ++) deg[i] = pstart[i+1]-pstart[i];

	ui *pend = new ui[n];

	for(ui i = 0;i < n;i ++) {
		ui &end = pend[i] = pstart[i];
		ui j = pstart[i+1];
		while(true) {
			while(end < j&&(deg[edges[end]] > deg[i]||(deg[edges[end]]==deg[i]&&edges[end]>i))) ++ end;
			while(j > end&&(deg[edges[j-1]] < deg[i]||(deg[edges[j-1]]==deg[i]&&edges[j-1]<i))) -- j;
			if(end >= j) break;
			swap(edges[end], edges[j-1]);
		}
		sort(edges+pstart[i], edges+pend[i]);
	}

	long long mtime1 = t.elapsed();

	ui *tri_cnt = new ui[m];
	memset(tri_cnt, 0, sizeof(ui)*m);

	memset(adj, 0, sizeof(ui)*n);

	long long res = 0;
	for(ui u = 0;u < n;u ++) {
		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = j+1;

		for(ui j = pstart[u];j < pend[u];j ++) {
			int v = edges[j];

			for(ui k = pstart[v];k < pend[v];k ++) {
				//++ tmp;
				if(adj[edges[k]]) {
					++ tri_cnt[j];
					++ tri_cnt[k];
					++ tri_cnt[adj[edges[k]]-1];
					res += 6;
				}
			}
		}

		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}

	long long mtime2 = t.elapsed();

	for(ui u = 0;u < n;u ++) {
		for(ui j = pend[u];j < pstart[u+1];j ++) {
			ui v = edges[j];

			ui r_id = binary_search(edges+pstart[v], pend[v]-pstart[v], u) + pstart[v];
			tri_cnt[j] = tri_cnt[r_id];
		}
	}

	long long mtime3 = t.elapsed();

	printf("#Triangles by TriC2: %lld\nTotal_Time: %lld\n", res, mtime3);
	printf("Sort_Time: %lld, CL_Time: %lld, Tri_Time: %lld\n", mtime1, mtime3-mtime2, mtime2-mtime1);

	delete[] tri_cnt;
	delete[] pend;
	delete[] adj;
}


void Graph::build_hashset() {
	if(hs != NULL) {
		printf("HashSet hs has already been initialized!\n");
		return ;
	}
#ifdef _HASHSET_
	hs = new HashSet*[n];
	//hs = new unordered_set<int>*[n];
#else
	hs = new hash_set<ui>*[n];
#endif

	for(ui i = 0;i < n;i ++) {
#ifdef _HASHSET_
		HashSet *lhs = hs[i] = new HashSet(pstart[i+1] - pstart[i]);
		//unordered_set<int> *lhs = hs[i] = new unordered_set<int>();
#else
		hash_set<ui> *lhs = hs[i] = new hash_set<ui>((pstart[i+1]-pstart[i])/2);
#ifdef _DENSE_
		lhs->set_empty_key(-1);
#endif
#endif
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(edges[j] > i) lhs->insert(edges[j]);
		}
	}
}

int Graph::find_hashset(ui a, ui b) {
#ifdef _HASHSET_
	if(a < b) return hs[a]->find(b) ;
	return hs[b]->find(a);
#else
	if(a < b) return hs[a]->find(b) != hs[a]->end();
	return hs[b]->find(a) != hs[b]->end();
#endif
}

void Graph::build_hashtable(ui *pend) {
#ifdef _HASHSET_
	printf("Wrong hashtable initialization!\n");
	return ;
#else
	if(ht != NULL) {
		printf("HashTable ht has already been initialized!\n");
		return ;
	}

	ht = new hash_map<ui, ui>*[n];

	for(ui i = 0;i < n;i ++) {
		ht[i] = new hash_map<ui, ui>(pend[i]-pstart[i]);
		hash_map<ui, ui> &lht = *ht[i];
#ifdef _DENSE_
		lht.set_empty_key(-1);
#endif
		for(ui j = pstart[i];j < pend[i];j ++) {
			lht[edges[j]] = 0;
		}
	}
#endif
}

void Graph::increment_hashtable(ui u, ui v) {
#ifdef _HASHSET_
	printf("Wrong hashtable search!\n");
	return ;
#else
	hash_map<ui, ui> &lht = *ht[u];
	++ lht[v];
#endif
}

ui Graph::find_hashtable(ui u, ui v) {
#ifdef _HASHSET_
	printf("Wrong hashtable search!\n");
	return 0;
#else
	hash_map<ui, ui> &lht = *ht[u];
	return lht[v];
#endif
}
