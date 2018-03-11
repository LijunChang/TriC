#include "HashSet.h"

HashSet::HashSet(int size) {
	P = size*LOAD;

	head = val = next = NULL;
	head = new ui[P];
	val = new ui[size];
	next = new ui[size];

	for(ui i = 0;i < P;i ++) head[i] = -1;

	capacity = size;
	cur = 0;
}

HashSet::~HashSet() {
	if(head != NULL) {
		delete[] head;
		head = NULL;
	}
	if(val != NULL) {
		delete[] val;
		val = NULL;
	}
	if(next != NULL) {
		delete[] next;
		next = NULL;
	}
}

void HashSet::insert(int v) {
#ifdef _DEBUG_
	if(cur == capacity) {
		printf("The hashset is already full.\n");
		return ;
	}
#endif

	int hv = v%P;
	val[cur] = v;
	next[cur] = head[hv];
	head[hv] = cur ++;
}

int HashSet::find(int v) {
	int hv = v%P;
	for(ui i = head[hv];i != -1;i = next[i]) {
		if(val[i] == v) return 1;
	}
	return 0;
}
