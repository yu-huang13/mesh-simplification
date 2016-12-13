//×îÐ¡¶Ñ
#ifndef HEAP_H
#define HEAP_H

#include <vector>
#include "pair.h"
using namespace std;

class Heap
{
public:
	Heap(vector<Pair>& arrayPairf, vector<int>& arrayPairHf);
	bool empty();
	int getMin();
	int delMin();
	void modify(int rank);
	void print();

private:
	inline bool inHeap(int i);
	inline int Parent(int i);
	inline int LChild(int i);
	inline int RChild(int i);
	inline int lastInternal();
	inline bool parentValid(int i);
	inline bool LChildValid(int i);
	inline bool RChildValid(int i);
	inline int Smaller(int i, int j);
	inline int properParent(int p);
	int percolateDown(int i);
	int percolateUp(int i);
	void heapify();

	vector<int>& Array;
	int size;

	vector<Pair>& pairs;
};

 


#endif