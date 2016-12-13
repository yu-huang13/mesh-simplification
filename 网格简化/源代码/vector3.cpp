#include <cmath>

#include "vector3.h"

const double EPS = 1e-6;

vector3::vector3(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}
vector3::vector3(const vector3 &v){
	x = v.x; y = v.y; z = v.z;
}

vector3 operator - (const vector3 &A){
	return vector3(-A.x, -A.y, -A.z);
}

vector3 operator + (const vector3 &A, const vector3 &B){
	return vector3(A.x + B.x, A.y + B.y, A.z + B.z);
}
vector3 operator - (const vector3 &A, const vector3 &B){
	return vector3(A.x - B.x, A.y - B.y, A.z - B.z);
}
vector3 operator * (const vector3 &A, const double &num){
	return vector3(A.x * num, A.y * num, A.z * num);
}
vector3 operator / (const vector3 &A, const double &num){
	return vector3(A.x / num, A.y / num, A.z / num);
}
vector3 operator * (const vector3 &A, const vector3 &B){
	return vector3(A.y * B.z - A.z * B.y, A.z * B.x - A.x * B.z, A.x * B.y - A.y * B.x);
}

vector3& operator += (vector3 &A, const vector3 &B){
	A = A + B;
	return A;
}
vector3& operator -= (vector3 &A, const vector3 &B){
	A = A - B;
	return A;
}
vector3& operator *= (vector3 &A, const double &num){
	A = A * num;
	return A;
}
vector3& operator /= (vector3 &A, const double &num){
	A = A / num;
	return A;
}
vector3& operator *= (vector3 &A, const vector3 &B){
	A = A * B;
	return A;
}

double vector3::dot(const vector3 &B){
	return x * B.x + y * B.y + z * B.z;
}
double vector3::module(){
	return sqrt(x * x + y * y + z * z);
}
double vector3::module2(){
	return x * x + y * y + z * z;
}
vector3 vector3::unitVector(){
	return (*this) / module();
}

bool vector3::isZero(){
	return fabs(x) < EPS && fabs(y) < EPS && fabs(z) < EPS;
}
vector3 vector3::reflect(const vector3 &N){
	return N * dot(N) * 2 - (*this);
}
bool vector3::refract(const vector3 &N, const double n, vector3& ref){//N为法向量
	vector3 V = unitVector();
	double cosH1 = V.dot(N);
	double temp = 1 - (1 - cosH1 * cosH1) / (n * n);
	if (temp > EPS){
		double cosH2 = sqrt(temp);
		ref = V * ((-1) / (n)) - N * (cosH2 - cosH1 / n);
		return true;
	}
	else
		return false;
}

bool operator == (const vector3 &A, const vector3 &B){
	return fabs(A.x - B.x) < EPS && fabs(A.y - B.y) < EPS && fabs(A.z - B.z) < EPS;
}

ostream& operator << (ostream& output, vector3 &A){
	output << setprecision(12) << A.x << " " << A.y << " " << A.z;
	return output;
}