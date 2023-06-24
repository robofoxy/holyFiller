#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <list>
#include <set>
#include <unordered_set>
#include <array>

using namespace std;

struct Vertex
{
	float* coords, * normals; //3d coordinates etc
	int idx; //who am i; verts[idx]

	vector< int > vertList; //adj vvertices;
	vector< int > triList;
	vector< int > edgeList;

	Vertex(int i, float* c) : idx(i), coords(c) {};
};

struct Edge
{
	int idx; //edges[idx]
	int v1i, v2i; //endpnts
	float length;
	Edge(int id, int v1, int v2) : idx(id), v1i(v1), v2i(v2) { computeLength(); };

	void computeLength()
	{
		length = 7;
	}
};

struct Triangle
{
	int idx; //tris[idx]
	int v1i, v2i, v3i;
	Triangle(int id, int v1, int v2, int v3) : idx(id), v1i(v1), v2i(v2), v3i(v3) {};
};

class Mesh
{
private:
	
	void addEdge(int v1, int v2);
	void addVertex(float x, float y, float z);
	bool makeVertsNeighbor(int v1i, int v2i);
public:
	vector< Vertex* > verts;
	vector<Eigen::Vector3f> transVerts;
	vector< Triangle* > tris;
	vector< Edge* > edges;
	void getYi(Eigen::MatrixXd&);
	void translateV(Eigen::MatrixXd);
	Mesh() {};
	void createCube(float side);
	void loadOff(char* name);
	void loadObj(const char* name);
	void addTriangle(int v1, int v2, int v3);
};

enum METHOD
{
	AREA,
	ANGLE
};

vector<int> holyFiller(Mesh*, enum METHOD method);
void boundaryLoopDetector(const std::vector<Triangle*>& tris, vector<vector<int>>& boundaryLoops);
void holyFillerHelper(Mesh*, const vector<int>& boundaryLoop, vector<int>& filled, enum METHOD method);

