#include "pair.h"

Pair::Pair(pair<int, int> E) : edge(E), deltaV(0), exist(true){}

size_t Hash_pair::operator () (const pair<int, int>& p) const {
	return p.first * 3 + p.second * 7 + 11;
}