#ifndef PAIR_H
#define PAIR_H

#include <utility>
#include <vector>
using namespace std;

struct Pair
{
	Pair(pair<int, int> E = pair<int, int>(0, 0));
	pair<int, int> edge;
	double deltaV;
	int rank_inHeap;
	bool exist;
};

struct Hash_pair{
	size_t operator () (const pair<int, int>& p) const;
};


#endif