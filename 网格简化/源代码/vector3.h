#ifndef VECTOR3_H
#define VECTOR3_H

#include <iostream>
#include <iomanip>

using namespace std;

struct vector3
{
	double x, y, z;
	vector3(double xx = 0, double yy = 0, double zz = 0);
	vector3(const vector3 &v);

	friend vector3 operator - (const vector3 &A);

	friend vector3 operator + (const vector3 &A, const vector3 &B);
	friend vector3 operator - (const vector3 &A, const vector3 &B);
	friend vector3 operator * (const vector3 &A, const double &num);
	friend vector3 operator / (const vector3 &A, const double &num);
	friend vector3 operator * (const vector3 &A, const vector3 &B);

	friend vector3& operator += (vector3 &A, const vector3 &B);
	friend vector3& operator -= (vector3 &A, const vector3 &B);
	friend vector3& operator *= (vector3 &A, const double &num);
	friend vector3& operator /= (vector3 &A, const double &num);
	friend vector3& operator *= (vector3 &A, const vector3 &B);

	friend bool operator == (const vector3 &A, const vector3 &B);

	friend ostream& operator << (ostream& output, vector3 &A);

	double dot(const vector3 &B);
	double module();
	double module2();
	vector3 unitVector();

	bool isZero();
	vector3 reflect(const vector3 &N);//N: unit normal vector3
	bool refract(const vector3 &N, const double n, vector3& ref);//N; unit normal vector3; n: relative index of refraction
};










#endif