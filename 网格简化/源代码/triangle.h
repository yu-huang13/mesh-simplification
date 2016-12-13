#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <iostream>
#include <iomanip>
using namespace std;

struct Triangle
{
	Triangle(int A, int B, int C);
	Triangle();

	friend ostream& operator << (ostream& output, Triangle& T);
	friend bool operator == (const Triangle& T1, const Triangle& T2);
	size_t operator () (const Triangle& T) const;//hashº¯Êý

	int p[3];
	bool exist;
};


#endif