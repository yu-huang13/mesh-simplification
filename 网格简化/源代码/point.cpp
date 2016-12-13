#include "point.h"

Point::Point(double x, double y, double z) : pos(x, y, z), exist(true) {}

Point::Point(vector3 v) : pos(v), exist(true) {}

Point::Point() : exist(true) {}

ostream& operator << (ostream& output, Point& P){
	output << P.pos;
	return output;
}