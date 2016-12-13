#include "triangle.h"


Triangle::Triangle(int A, int B, int C){
	p[0] = A; p[1] = B; p[2] = C; exist = true;
}

Triangle::Triangle() : exist(true) {}

ostream& operator << (ostream& output, Triangle& T){
	output << T.p[0] << " " << T.p[1] << " " << T.p[2];
	return output;
}

bool operator == (const Triangle& T1, const Triangle& T2){
	return T1.p[0] == T2.p[0] && T1.p[1] == T2.p[1] && T1.p[2] == T2.p[2];
}

size_t Triangle::operator () (const Triangle& T) const{//hashº¯Êý
	return p[0] * 2 + p[1] * 3 + p[2] * 5 + 11;
}