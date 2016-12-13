#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <unordered_set>

#include "point.h"
#include "triangle.h"
#include "vector3.h"
#include "heap.h"

using namespace std;

const double EPS = 1e-200;

vector<Point> points; // 1 - (size - 1)  所有点的坐标
double*** pointQ;//1 - (size - 1) pointQ[i]:points[i]对应的矩阵
vector<Triangle> triangles; // 1 - (size - 1)  所有三角形
vector < vector<int> > pointTriangle;//1 - (size - 1)  [i]: points[i]关联的三角形
vector < vector<int> > pointPairs; //1 - (size - 1)  [i]:points[i]关联的pair
vector<Pair> pairs;//0 - (size - 1)  所有边
vector<int> pairsH;//0 - (size - 1)  堆

double* temp_4Array;
double** temp_44Array;

int all_triangle;
int triangleNum;

//将a、b从小到大排序
inline void SortPair(int& a, int& b){
	if (a > b) swap(a, b);
}

//将a、b、c从小到大排序
void SortTri(int &a, int &b, int&c){
	if (c < a) swap(c, a);
	if (c < b) swap(c, b);
	if (b < a) swap(b, a);
}

//new出4*4二维数组
void new_44Array(double** &a){
	a = new double*[4];
	for (int i = 0; i < 4; ++i)
		a[i] = new double[4];
}

//delete4*4二维数组
void delete_44Array(double** &a){
	for (int i = 0; i < 4; ++i)
		delete[]a[i];
	delete[]a;
}

//delete三维数组
void delete_3d_44Array(double*** &a, int size){
	for (int i = 0; i < size; ++i)
		delete_44Array(a[i]);
	delete[]a;
}

//K1 = K1 + K2
void add_44Array(double** K1, double** K2){
	for (int i = 0; i < 4; ++i)
	for (int j = 0; j < 4; ++j)
		K1[i][j] = K1[i][j] + K2[i][j];
}

//result = K1 + K2
void add_44Array(double** K1, double** K2, double** result){
	for (int i = 0; i < 4; ++i)
	for (int j = 0; j < 4; ++j)
		result[i][j] = K1[i][j] + K2[i][j];
}

//将4*4矩阵初始化为0
void init0_44Array(double** a){
	for (int i = 0; i < 4; ++i)
	for (int j = 0; j < 4; ++j)
		a[i][j] = 0;
}

//读取obj文件
bool readFile(char* filename)
{
	ifstream fin(filename);
	if (!fin){
		cout << "fail to open " << filename << endl;
		return false;
	}
	string line;
	string mark;
	stringstream fin2;
	double P[3];
	int T[3];
	points.push_back(Point(0, 0, 0));//为使points的下标从1开始有效，预先加入
	triangles.push_back(Triangle());//为使Triangle的下标从1开始有效，预先加入
	while (getline(fin, line)){
		fin2.clear();
		fin2.str(line);
		fin2 >> mark;
		if (mark == "v"){
			for (int i = 0; i < 3; ++i)
				fin2 >> P[i];
			points.push_back(Point(P[0], P[1], P[2]));
		}
		if (mark == "f"){
			for (int i = 0; i < 3; ++i)
				fin2 >> T[i];
			SortTri(T[0], T[1], T[2]);
			triangles.push_back(Triangle(T[0], T[1], T[2]));
		}
	}
	all_triangle = triangleNum = triangles.size() - 1;

	fin.close();
	return true;
}

//初始化pointTriangle数组
void init_pointTriangle()
{
	pointTriangle.resize(points.size());
	int trisize = triangles.size();
	for (int i = 1; i < trisize; ++i){
		for (int k = 0; k < 3; ++k)
			pointTriangle[triangles[i].p[k]].push_back(i);
	}
}

//初始化pairs数组
void init_Pair()
{
	pointPairs.resize(points.size());
	int trisize = triangles.size();
	int a, b;//a < b
	pair<int, int> temp;
	unordered_set<pair<int, int>, Hash_pair> recordP;
	for (int i = 1; i < trisize; ++i){
		for (int k = 0; k < 3; ++k){
			a = triangles[i].p[k]; b = triangles[i].p[(k + 1) % 3]; if (a > b) swap(a, b);
			temp = pair<int, int>(a, b);
			if (recordP.find(temp) == recordP.end()){//没找到
				recordP.emplace(temp);
				pointPairs[a].push_back(pairs.size());
				pointPairs[b].push_back(pairs.size());
				pairs.push_back(Pair(temp));	
			}
		}
	}
	int pairsize = pairs.size();
	for (int i = 0; i < pairsize; ++i){
		pairsH.push_back(i);
		pairs[i].rank_inHeap = i;
	}
}

//计算三角形ABC的平面方程四个参数，存在p中
void cal_abcd(const vector3& A, const vector3& B, const vector3& C, double *p){
	vector3 u1 = B - A;
	vector3 u2 = C - A;
	vector3 N = (u1 * u2).unitVector();
	p[0] = N.x; p[1] = N.y; p[2] = N.z;
	p[3] = -N.dot(A);
}

//Kp = p1T * p1
void cal_Kp(double** Kp, double* p1, double* p2){//p1的转置乘p2
	for (int i = 0; i < 4; ++i)
	for (int j = 0; j < 4; ++j)
		Kp[i][j] = p1[i] * p2[j];
}

//计算point[pointRank]的误差矩阵，结果存在Q中
void cal_Q(double** Q, int pointRank){
	init0_44Array(Q);
	double *p = new double[4];
	for (vector<int>::iterator it = pointTriangle[pointRank].begin(); it != pointTriangle[pointRank].end(); ++it){
		if (!triangles[*it].exist) continue;
		cal_abcd(points[triangles[*it].p[0]].pos, points[triangles[*it].p[1]].pos, points[triangles[*it].p[2]].pos, p);
		cal_Kp(temp_44Array, p, p);
		add_44Array(Q, temp_44Array);
	}
	delete[]p;
}

//初始化pointQ数组
void init_pointQ(){
	int size = points.size();
	pointQ = new double**[size];
	for (int i = 0; i < size; ++i){
		new_44Array(pointQ[i]);
	}
	for (int i = 1; i < size; ++i){//计算point[i]的误差矩阵
		cal_Q(pointQ[i], i);
	}
}

//v * Q * vT
double cal_VQV(double* v, double** Q){
	double result = 0;
	double temp;
	for (int i = 0; i < 4; ++i){
		temp = 0;
		for (int j = 0; j < 4; ++j){
			temp += Q[i][j] * v[j];
		}
		result += v[i] * temp;
	}
	return result;
}

//计算3 * 3矩阵行列式
double det33(const vector3& v1, const vector3& v2, const vector3& v3){//列向量
	return v1.x * (v2.y * v3.z - v3.y * v2.z) - v2.x * (v1.y * v3.z - v3.y * v1.z) + v3.x * (v1.y * v2.z - v2.y * v1.z);
}

//计算收缩点（直接取中点）
vector3 getShrinkPoint(int A, int B){
	return (points[A].pos + points[B].pos) / 2;
}

//计算收缩点（误差矩阵）
vector3 getShrinkPoint(int A, int B, double** Q){
	add_44Array(pointQ[A], pointQ[B], Q);
	vector3 v1, v2, v3, v4;
	v1.x = Q[0][0]; v2.x = Q[0][1]; v3.x = Q[0][2]; v4.x = -Q[0][3];
	v1.y = Q[1][0]; v2.y = Q[1][1]; v3.y = Q[1][2]; v4.y = -Q[1][3];
	v1.z = Q[2][0]; v2.z = Q[2][1]; v3.z = Q[2][2]; v4.z = -Q[2][3];
	double det = det33(v1, v2, v3);
	if (fabs(det) < EPS) {
		return getShrinkPoint(A, B);
	}
	return vector3(det33(v4, v2, v3) / det, det33(v1, v4, v3) / det, det33(v1, v2, v4) / det);
}

//计算pair的deltaV
double cal_deltaV(const pair<int, int>& edge){
	add_44Array(pointQ[edge.first], pointQ[edge.second], temp_44Array);
	vector3 shrinkP = getShrinkPoint(edge.first, edge.second, temp_44Array);
	temp_4Array[0] = shrinkP.x; temp_4Array[1] = shrinkP.y; temp_4Array[2] = shrinkP.z; temp_4Array[3] = 1;
	return cal_VQV(temp_4Array, temp_44Array);
}

//初始化每条边的deltaV
void init_deltaV(){
	int size = pairs.size();
	vector3 temp;
	for (int i = 0; i < size; ++i){
		pairs[i].deltaV = cal_deltaV(pairs[i].edge);
	}
}

//判断triangles[triRank]是否有顶点pRank
bool triangle_has_point(int triRank, int pRank){
	for (int i = 0; i < 3; ++i){
		if (triangles[triRank].p[i] == pRank)
			return true;
	}
	return false;
}

//缩点, v、temp为暂存量,以避免反复new与delete
void Shrink(Heap& heap)
{
	int minRank;
	
	while (!heap.empty()){//出堆
		minRank = heap.delMin();
		if (pairs[minRank].exist){
			//cout<<"mimRank: " << pairs[minRank].deltaV << endl;
			break;
		}
	}
	if (heap.empty()){ cout << "heap empty()" << endl; return; }

	int x = pairs[minRank].edge.first, y = pairs[minRank].edge.second;//将y收缩到x
	points[y].exist = false;

	//确定收缩点的位置
	//points[x].pos = getShrinkPoint(x, y);
	points[x].pos = getShrinkPoint(x, y, temp_44Array);

	//更新三角形
	for (vector<int>::iterator it = pointTriangle[y].begin(); it != pointTriangle[y].end(); ++it){
		if (!triangles[*it].exist) continue;
		if (triangle_has_point(*it, x)){
			triangles[*it].exist = false;
			--triangleNum;
		}
		else{
			triangles[*it].p[0] == y ? triangles[*it].p[0] = x :
				triangles[*it].p[1] == y ? triangles[*it].p[1] = x : triangles[*it].p[2] = x;
			SortTri(triangles[*it].p[0], triangles[*it].p[1], triangles[*it].p[2]);
			pointTriangle[x].push_back(*it);			
		}
	}

	//三角形判重
	unordered_set<Triangle, Triangle> recordT;
	for (vector<int>::iterator it = pointTriangle[x].begin(); it != pointTriangle[x].end();++it){
		if (!triangles[*it].exist) continue;
		if (recordT.find(triangles[*it]) == recordT.end())//没找到
			recordT.emplace(triangles[*it]);
		else{//找到
			triangles[*it].exist = false;
			--triangleNum;
		}
	}

	//更新pairs
	for (vector<int>::iterator it = pointPairs[y].begin(); it != pointPairs[y].end(); ++it){
		if (!pairs[*it].exist) continue;
		pairs[*it].edge.first == y ? pairs[*it].edge.first = x : pairs[*it].edge.second = x;
		if (pairs[*it].edge.first == pairs[*it].edge.second){
			pairs[*it].exist = false;
		}
		else{
			SortPair(pairs[*it].edge.first, pairs[*it].edge.second);
			pointPairs[x].push_back(*it);
		}
			
	}

	//pair判重
	unordered_set<pair<int, int>, Hash_pair> recordP;
	for (vector<int>::iterator it = pointPairs[x].begin(); it != pointPairs[x].end(); ++it){
		if (!pairs[*it].exist) continue;
		if (recordP.find(pairs[*it].edge) == recordP.end())//没找到
			recordP.emplace (pairs[*it].edge);
		else//找到
			pairs[*it].exist = false;
	}

	//确定需要更新Q的点
	for (vector<int>::iterator it = pointPairs[x].begin(); it != pointPairs[x].end(); ++it){
		if (!pairs[*it].exist) continue;
		if (pairs[*it].edge.first == x){
			cal_Q(pointQ[pairs[*it].edge.second], pairs[*it].edge.second);
		}
		else{
			cal_Q(pointQ[pairs[*it].edge.first], pairs[*it].edge.first);
		}	
	}
	cal_Q(pointQ[x], x);


	//更新pair的deltaV
	//优化：hash判重
	unordered_set<pair<int, int>, Hash_pair> updateP;
	for (vector<int>::iterator it = pointPairs[x].begin(); it != pointPairs[x].end(); ++it){
		if (!pairs[*it].exist) continue;
		int rankP = pairs[*it].edge.first == x ? pairs[*it].edge.second : pairs[*it].edge.first;
		for (vector<int>::iterator itU = pointPairs[rankP].begin(); itU != pointPairs[rankP].end(); ++itU){
			if (!pairs[*itU].exist) continue;
			if (updateP.find(pairs[*itU].edge) == updateP.end()){//没找到
				updateP.insert(pairs[*itU].edge);
				pairs[*itU].deltaV = cal_deltaV(pairs[*itU].edge);
				heap.modify(pairs[*itU].rank_inHeap);
			}
		}
	}
}

//网格简化，rate为百分比
void MeshSimplification(Heap& heap, double rate){
	int target_num = all_triangle * rate;
	while (triangleNum > target_num){
		Shrink(heap);
	}
}

//输出简化后的模型
void print_obj(char* filename){
	ofstream fout(filename);
	int size = points.size();
	int* mark = new int[points.size()];
	int rank = 1;
	for (int i = 1; i < size; ++i){
		if (!points[i].exist) continue;
		mark[i] = rank;
		++rank;
		fout <<"v "<< points[i].pos<<endl;
	}
	size = triangles.size();
	for (int i = 1; i < size; ++i){
		if (triangles[i].exist)
			fout <<"f "<< mark[triangles[i].p[0]] << " " << mark[triangles[i].p[1]] << " " << mark[triangles[i].p[2]] << endl;
	}
	delete[]mark;
	fout.close();
}

//格式：main.exe 输入文件名 输出文件名 简化百分比（三角形面片数量）
int main(int argc, char** argv)
{
	temp_4Array = new double[4];
	new_44Array(temp_44Array);

	readFile(argv[1]);
	init_pointTriangle();
	init_Pair();
	init_pointQ();
	init_deltaV();
	Heap heap(pairs, pairsH);
	MeshSimplification(heap, atof(argv[3]));
	print_obj(argv[2]);
	
	delete_3d_44Array(pointQ, points.size());
	delete[]temp_4Array;
	delete_44Array(temp_44Array);

	system("pause");
	return 0;
}
