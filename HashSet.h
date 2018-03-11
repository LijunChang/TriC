#ifndef _HASH_SET_
#define _HASH_SET_

#include "Utility.h"

class HashSet {
private:
	ui *head;
	ui *val;
	ui *next;
	int cur, capacity, P;

public:
	HashSet(int size) ;
	~HashSet() ;

	void insert(int v) ;
	int find(int v) ;
};

#endif
