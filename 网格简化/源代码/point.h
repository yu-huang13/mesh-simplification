#ifndef POINT_H
#define POINT_H

#include <iostream>
#include "vector3.h"
using namespace std;

struct Point
{
	Point(double x, double y, double z);
	Point(vector3 v);
	Point();

	friend ostream& operator << (ostream& output, Point& P);

	vector3 pos;
	bool exist;
};



#endif