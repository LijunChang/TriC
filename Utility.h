#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <queue>
#include <set>

#define NDEBUG
#include <cassert>

#define _HASHSET_

#ifndef _HASHSET_
#define _DENSE_
#ifdef _DENSE_
	#define hash_set dense_hash_set
	#define hash_map dense_hash_map
#else
	#define hash_set sparse_hash_set
	#define hash_map sparse_hash_map
#endif

	#include <google/dense_hash_set>
	#include <google/dense_hash_map>
	using google::dense_hash_set;
	using google::dense_hash_map;

	#include <google/sparse_hash_set>
	#include <google/sparse_hash_map>
	using google::sparse_hash_set;
	using google::sparse_hash_map;
#endif

typedef unsigned int ui;
#define pb push_back
#define mp make_pair

const int LOAD = 2;

FILE *open_file(const char *file_name, const char *mode) ;

#endif
