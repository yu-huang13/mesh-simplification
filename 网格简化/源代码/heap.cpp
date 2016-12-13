#include "heap.h"
#include <iostream>
using namespace std;
Heap::Heap(vector<Pair>& pairsf, vector<int>& Arrayf) : pairs(pairsf), Array(Arrayf){
	size = Arrayf.size();
	heapify();
}


bool Heap::empty(){
	return size == 0;
}


inline bool Heap::inHeap(int i){
	return i > -1 && i < size;
}


inline int Heap::Parent(int i){
	return (i - 1) / 2;
}


inline int Heap::LChild(int i){
	return 2 * i + 1;
}


inline int Heap::RChild(int i){
	return 2 * i + 2;
}

inline int Heap::lastInternal(){
	return Parent(size - 1);
}


inline bool Heap::parentValid(int i){
	return i > 0;
}


inline bool Heap::LChildValid(int i){
	return inHeap(LChild(i));
}


inline bool Heap::RChildValid(int i){
	return inHeap(RChild(i));
}


inline int Heap::Smaller(int i, int j){//等时前者优先
	return pairs[Array[j]].deltaV < pairs[Array[i]].deltaV ? j : i;
}


inline int Heap::properParent(int p){
	return RChildValid(p) ? Smaller(Smaller(p, LChild(p)), RChild(p)) :
		(LChildValid(p) ? Smaller(p, LChild(p)) : p);
}


int Heap::percolateDown(int i){
	int j;
	while (i != (j = properParent(i))){
		swap(pairs[Array[i]].rank_inHeap, pairs[Array[j]].rank_inHeap);
		swap(Array[i], Array[j]);
		i = j;
	}
	return i;
}


int Heap::percolateUp(int i){
	while (parentValid(i)){
		int j = Parent(i);
		if (j == Smaller(i, j)) break;
		swap(pairs[Array[i]].rank_inHeap, pairs[Array[j]].rank_inHeap);
		swap(Array[i], Array[j]);
		i = j;
	}
	return i;
}

int Heap::getMin(){//调用前应保证堆不为空
	return Array[0];
}


int Heap::delMin(){//调用前应保证堆不为空
	int minT = Array[0];
	Array[0] = Array[--size]; Array.pop_back();
	percolateDown(0);
	return minT;
}

void Heap::heapify(){
	for (int i = lastInternal(); inHeap(i); --i){
		percolateDown(i);
	}
}

void Heap::modify(int rank){
	percolateUp(rank);
	percolateDown(rank);
}

void Heap::print(){
	int j = 2;
	cout << pairs[Array[0]].deltaV << endl;
	while (true){
		for (int i = j - 1; i < 2 * j - 1; ++i){
			if (i >= size) return;
			cout << pairs[Array[i]].deltaV << " ";
		}
		j = j * 2;
		cout << endl;
	}
	
}