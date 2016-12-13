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

vector<Point> points; // 1 - (size - 1)  ���е������
double*** pointQ;//1 - (size - 1) pointQ[i]:points[i]��Ӧ�ľ���
vector<Triangle> triangles; // 1 - (size - 1)  ����������
vector < vector<int> > pointTriangle;//1 - (size - 1)  [i]: points[i]������������
vector < vector<int> > pointPairs; //1 - (size - 1)  [i]:points[i]������pair
vector<Pair> pairs;//0 - (size - 1)  ���б�
vector<int> pairsH;//0 - (size - 1)  ��

double* temp_4Array;
double** temp_44Array;

int all_triangle;
int triangleNum;

//��a��b��С��������
inline void SortPair(int& a, int& b){
	if (a > b) swap(a, b);
}

//��a��b��c��С��������
void SortTri(int &a, int &b, int&c){
	if (c < a) swap(c, a);
	if (c < b) swap(c, b);
	if (b < a) swap(b, a);
}

//new��4*4��ά����
void new_44Array(double** &a){
	a = new double*[4];
	for (int i = 0; i < 4; ++i)
		a[i] = new double[4];
}

//delete4*4��ά����
void delete_44Array(double** &a){
	for (int i = 0; i < 4; ++i)
		delete[]a[i];
	delete[]a;
}

//delete��ά����
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

//��4*4�����ʼ��Ϊ0
void init0_44Array(double** a){
	for (int i = 0; i < 4; ++i)
	for (int j = 0; j < 4; ++j)
		a[i][j] = 0;
}

//��ȡobj�ļ�
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
	points.push_back(Point(0, 0, 0));//Ϊʹpoints���±��1��ʼ��Ч��Ԥ�ȼ���
	triangles.push_back(Triangle());//ΪʹTriangle���±��1��ʼ��Ч��Ԥ�ȼ���
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

//��ʼ��pointTriangle����
void init_pointTriangle()
{
	pointTriangle.resize(points.size());
	int trisize = triangles.size();
	for (int i = 1; i < trisize; ++i){
		for (int k = 0; k < 3; ++k)
			pointTriangle[triangles[i].p[k]].push_back(i);
	}
}

//��ʼ��pairs����
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
			if (recordP.find(temp) == recordP.end()){//û�ҵ�
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

//����������ABC��ƽ�淽���ĸ�����������p��
void cal_abcd(const vector3& A, const vector3& B, const vector3& C, double *p){
	vector3 u1 = B - A;
	vector3 u2 = C - A;
	vector3 N = (u1 * u2).unitVector();
	p[0] = N.x; p[1] = N.y; p[2] = N.z;
	p[3] = -N.dot(A);
}

//Kp = p1T * p1
void cal_Kp(double** Kp, double* p1, double* p2){//p1��ת�ó�p2
	for (int i = 0; i < 4; ++i)
	for (int j = 0; j < 4; ++j)
		Kp[i][j] = p1[i] * p2[j];
}

//����point[pointRank]�������󣬽������Q��
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

//��ʼ��pointQ����
void init_pointQ(){
	int size = points.size();
	pointQ = new double**[size];
	for (int i = 0; i < size; ++i){
		new_44Array(pointQ[i]);
	}
	for (int i = 1; i < size; ++i){//����point[i]��������
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

//����3 * 3��������ʽ
double det33(const vector3& v1, const vector3& v2, const vector3& v3){//������
	return v1.x * (v2.y * v3.z - v3.y * v2.z) - v2.x * (v1.y * v3.z - v3.y * v1.z) + v3.x * (v1.y * v2.z - v2.y * v1.z);
}

//���������㣨ֱ��ȡ�е㣩
vector3 getShrinkPoint(int A, int B){
	return (points[A].pos + points[B].pos) / 2;
}

//���������㣨������
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

//����pair��deltaV
double cal_deltaV(const pair<int, int>& edge){
	add_44Array(pointQ[edge.first], pointQ[edge.second], temp_44Array);
	vector3 shrinkP = getShrinkPoint(edge.first, edge.second, temp_44Array);
	temp_4Array[0] = shrinkP.x; temp_4Array[1] = shrinkP.y; temp_4Array[2] = shrinkP.z; temp_4Array[3] = 1;
	return cal_VQV(temp_4Array, temp_44Array);
}

//��ʼ��ÿ���ߵ�deltaV
void init_deltaV(){
	int size = pairs.size();
	vector3 temp;
	for (int i = 0; i < size; ++i){
		pairs[i].deltaV = cal_deltaV(pairs[i].edge);
	}
}

//�ж�triangles[triRank]�Ƿ��ж���pRank
bool triangle_has_point(int triRank, int pRank){
	for (int i = 0; i < 3; ++i){
		if (triangles[triRank].p[i] == pRank)
			return true;
	}
	return false;
}

//����, v��tempΪ�ݴ���,�Ա��ⷴ��new��delete
void Shrink(Heap& heap)
{
	int minRank;
	
	while (!heap.empty()){//����
		minRank = heap.delMin();
		if (pairs[minRank].exist){
			//cout<<"mimRank: " << pairs[minRank].deltaV << endl;
			break;
		}
	}
	if (heap.empty()){ cout << "heap empty()" << endl; return; }

	int x = pairs[minRank].edge.first, y = pairs[minRank].edge.second;//��y������x
	points[y].exist = false;

	//ȷ���������λ��
	//points[x].pos = getShrinkPoint(x, y);
	points[x].pos = getShrinkPoint(x, y, temp_44Array);

	//����������
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

	//����������
	unordered_set<Triangle, Triangle> recordT;
	for (vector<int>::iterator it = pointTriangle[x].begin(); it != pointTriangle[x].end();++it){
		if (!triangles[*it].exist) continue;
		if (recordT.find(triangles[*it]) == recordT.end())//û�ҵ�
			recordT.emplace(triangles[*it]);
		else{//�ҵ�
			triangles[*it].exist = false;
			--triangleNum;
		}
	}

	//����pairs
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

	//pair����
	unordered_set<pair<int, int>, Hash_pair> recordP;
	for (vector<int>::iterator it = pointPairs[x].begin(); it != pointPairs[x].end(); ++it){
		if (!pairs[*it].exist) continue;
		if (recordP.find(pairs[*it].edge) == recordP.end())//û�ҵ�
			recordP.emplace (pairs[*it].edge);
		else//�ҵ�
			pairs[*it].exist = false;
	}

	//ȷ����Ҫ����Q�ĵ�
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


	//����pair��deltaV
	//�Ż���hash����
	unordered_set<pair<int, int>, Hash_pair> updateP;
	for (vector<int>::iterator it = pointPairs[x].begin(); it != pointPairs[x].end(); ++it){
		if (!pairs[*it].exist) continue;
		int rankP = pairs[*it].edge.first == x ? pairs[*it].edge.second : pairs[*it].edge.first;
		for (vector<int>::iterator itU = pointPairs[rankP].begin(); itU != pointPairs[rankP].end(); ++itU){
			if (!pairs[*itU].exist) continue;
			if (updateP.find(pairs[*itU].edge) == updateP.end()){//û�ҵ�
				updateP.insert(pairs[*itU].edge);
				pairs[*itU].deltaV = cal_deltaV(pairs[*itU].edge);
				heap.modify(pairs[*itU].rank_inHeap);
			}
		}
	}
}

//����򻯣�rateΪ�ٷֱ�
void MeshSimplification(Heap& heap, double rate){
	int target_num = all_triangle * rate;
	while (triangleNum > target_num){
		Shrink(heap);
	}
}

//����򻯺��ģ��
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

//��ʽ��main.exe �����ļ��� ����ļ��� �򻯰ٷֱȣ���������Ƭ������
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
