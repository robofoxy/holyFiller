#include "Mesh.h"

#define ALPHA sqrt(2)
#define INFINITE_PREVENTER 20

bool circumSphereCheck(Mesh* mesh, int v1_, int v2_, int v3_, int op_);
bool swapCheck(const vector<Eigen::Vector3i>& triangles, int iv1, int iv2, int op2, int op1);

struct edge {
	int X;
	int Y;

	edge() : X(0), Y(0) {};
	edge(const int& x, const int& y) : X(x), Y(y) {};
	edge(const edge& other) {
		X = other.X;
		Y = other.Y;
	};

	edge& operator=(const edge& other) {
		X = other.X;
		Y = other.Y;
		return *this;
	};

	bool operator==(const edge& other) const {
		if (X == other.X && Y == other.Y)
			return true;
		return false;
	};

	bool operator<(const edge& other) {
		if (X < other.X)
			return true;
		else if (X == other.X && Y == other.Y)
			return true;

		return false;
	};

	// this could be moved in to std::hash<edge>::operator()
	size_t operator()(const edge& pointToHash) const noexcept {
		size_t hash = pointToHash.X + 10 * pointToHash.Y;
		return hash;
	};

};

namespace std {
	template<> struct hash<edge>
	{
		std::size_t operator()(const edge& p) const noexcept
		{
			return p(p);
		}
	};
}


void boundaryLoopDetector(const vector<Triangle*>& tris, vector<vector<int>>& boundaryLoops)
{
	unordered_set<edge> edges;

	for (int i = 0; i < tris.size(); i++)
	{
		edge ed { tris[i]->v2i, tris[i]->v1i };
		unordered_set<edge>::iterator aEd = edges.find(ed);
		
		if (aEd == edges.end())
		{
			edge e{ tris[i]->v1i, tris[i]->v2i };
			edges.insert(e);
		}
		else
		{
			edges.erase(aEd);
		}

		edge ed1{ tris[i]->v3i, tris[i]->v2i };

		if((aEd = edges.find(ed1)) == edges.end())
		{
			edge e{ tris[i]->v2i, tris[i]->v3i };
			edges.insert(e);
		}
		else
		{
			edges.erase(aEd);
		}

		edge ed2{ tris[i]->v1i, tris[i]->v3i };

		if((aEd = edges.find(ed2)) == edges.end())
		{
			edge e{ tris[i]->v3i, tris[i]->v1i };
			edges.insert(e);
		}
		else
		{
			edges.erase(aEd);
		}
	}
	unordered_map<int, int> boundary;

	for (unordered_set<edge>::iterator ed = edges.begin(); ed != edges.end(); ed++)
	{
		boundary[(*ed).X] = (*ed).Y;
	}

	vector<int> boundaryLoop(boundary.size(), 0);
	
	int vi = -1;
	
	while (true)
	{
		if (vi == -1)
		{
			if (boundary.empty())
			{
				break;
			}
			
			vi = boundary.begin()->first;
			boundaryLoop.clear();
			boundaryLoop.push_back(vi);
		}
		else
		{
			int nVi = boundary[vi];
			boundary.erase(vi);
			
			if (nVi == boundaryLoop[0])
			{
				reverse(boundaryLoop.begin(), boundaryLoop.end());
				boundaryLoops.push_back(boundaryLoop);
				vi = -1;
			}
			else
			{
				boundaryLoop.push_back(nVi);
				vi = nVi;
			}
		}
	}
}

vector<Eigen::Vector3i> trianglesToBeInserted(Mesh* mesh, const vector<int>& boundaryLoop, const vector<vector<int>>& lambdas, vector<Eigen::Vector3d>& centroids)
{
	vector<pair<int, int>> sections;

	sections.push_back(pair<int, int>{ 0, boundaryLoop.size() - 1 });

	vector<Eigen::Vector3i> to_be;

	while (sections.size() > 0)
	{
		pair<int, int> section = sections.back();
		sections.pop_back();
		int k;

		if (section.second - section.first == 2)
		{
			k = section.first + 1;
		}
		else
		{
			k = lambdas[section.second - section.first - 1][section.first];
		}

		Eigen::Vector3i tr(boundaryLoop[section.first], boundaryLoop[k], boundaryLoop[section.second]);

		to_be.push_back(tr);

		Eigen::Vector3d vc;
		vc << (mesh->verts[boundaryLoop[section.first]]->coords[0] + mesh->verts[boundaryLoop[k]]->coords[0] + mesh->verts[boundaryLoop[section.second]]->coords[0]) / 3.0,
			(mesh->verts[boundaryLoop[section.first]]->coords[1] + mesh->verts[boundaryLoop[k]]->coords[1] + mesh->verts[boundaryLoop[section.second]]->coords[1]) / 3.0,
			(mesh->verts[boundaryLoop[section.first]]->coords[2] + mesh->verts[boundaryLoop[k]]->coords[2] + mesh->verts[boundaryLoop[section.second]]->coords[2]) / 3.0;
		centroids.push_back(vc);

		if (k - section.first > 1)
		{
			sections.push_back(pair<int, int>{ section.first, k });
		}

		if (section.second - k > 1)
		{
			sections.push_back(pair<int, int>{ k, section.second });
		}
	}

	return to_be;
}


void insertTriangles(Mesh* mesh, const vector<Eigen::Vector3i>& triangles)
{
	for (int i = 0; i < triangles.size(); i++)
	{
		mesh->addTriangle(triangles[i](0), triangles[i](1), triangles[i](2));
	}
}

vector<int> findOrigins(vector<int> face, int size) {
	sort(face.begin(), face.end());
	
	if (face[0] == INT_MIN)
	{
		if (face[1] == 0 && face[2] == size - 1)
		{
			return vector<int>{ size - 1, INT_MIN };
		}
		
		if (face[1] + 1 == face[2])
		{
			return vector<int>{ face[1], INT_MIN };
		}

		return vector<int>{ INT_MIN, INT_MIN };
	}

	if (face[0] == 0 && face[2] == size - 1)
	{
		if (face[1] == 1)
		{
			return vector<int>{ size - 1, 0 };
		}
		if (face[1] == size - 2)
		{
			return vector<int>{ size - 2, size - 1 };
		}
	}

	return vector<int>{ face[0], face[1] };
}

double trArea(const vector<Vertex*>& vertices, int i, int j, int k)
{
	Eigen::Vector3d a;
	a << vertices[j]->coords[0] - vertices[i]->coords[0], vertices[j]->coords[1] - vertices[i]->coords[1], vertices[j]->coords[2] - vertices[i]->coords[2];
	Eigen::Vector3d b;
	b << vertices[k]->coords[0] - vertices[i]->coords[0], vertices[k]->coords[1] - vertices[i]->coords[1], vertices[k]->coords[2] - vertices[i]->coords[2];

	return a.cross(b).norm() / 2.0;
}

Eigen::Vector3d trNormal(const vector<Vertex*>& vertices, int i, int j, int k)
{
	Eigen::Vector3d a;
	a << vertices[j]->coords[0] - vertices[i]->coords[0], vertices[j]->coords[1] - vertices[i]->coords[1], vertices[j]->coords[2] - vertices[i]->coords[2];
	Eigen::Vector3d b;
	b << vertices[k]->coords[0] - vertices[i]->coords[0], vertices[k]->coords[1] - vertices[i]->coords[1], vertices[k]->coords[2] - vertices[i]->coords[2];

	return a.cross(b).normalized();
}

void findEdgeTriangleNormals(Mesh* mesh, vector<vector<Eigen::Vector3d>>& edgeTriangleNormals, const vector<int>& boundaryLoop)
{
	unordered_set<int> boundaryVertices;
	
	for (int i = 0; i < boundaryLoop.size(); i++)
	{
		boundaryVertices.insert(boundaryLoop[i]);
	}

	for (int i = boundaryLoop.size() - 1; i > 0; i--)
	{
		if (i < boundaryLoop.size() - 1)
		{
			vector<Eigen::Vector3d> vec(i);
			edgeTriangleNormals.push_back(vec);
		}
		else
		{
			vector<Eigen::Vector3d> vec(boundaryLoop.size());
			edgeTriangleNormals.push_back(vec);
		}
	}

	for (int i = 0; i < mesh->tris.size(); i++)
	{
		int count = 0;
		vector<int> face(3, INT_MIN);
		unordered_set<int>::iterator it;
		if ((it = boundaryVertices.find(mesh->tris[i]->v1i)) != boundaryVertices.end()) { count++; face[0] = find(boundaryLoop.begin(), boundaryLoop.end(), *it) - boundaryLoop.begin(); }
		if ((it = boundaryVertices.find(mesh->tris[i]->v2i)) != boundaryVertices.end()) { count++; face[1] = find(boundaryLoop.begin(), boundaryLoop.end(), *it) - boundaryLoop.begin(); }
		if ((it = boundaryVertices.find(mesh->tris[i]->v3i)) != boundaryVertices.end()) { count++; face[2] = find(boundaryLoop.begin(), boundaryLoop.end(), *it) - boundaryLoop.begin(); }

		if (count > 1)
		{
			Eigen::Vector3d normal = trNormal(mesh->verts, mesh->tris[i]->v1i, mesh->tris[i]->v2i, mesh->tris[i]->v3i);
			vector<int> origins = findOrigins(face, boundaryLoop.size());

			if (origins[0] != INT_MIN)
			{
				edgeTriangleNormals[0][origins[0]] = normal;
			}

			if (origins[1] != INT_MIN)
			{
				edgeTriangleNormals[0][origins[1]] = normal;
			}
		}
	}
	for (int i = 0; i < boundaryLoop.size() - 2; i++)
	{
		edgeTriangleNormals[1][i] = trNormal(mesh->verts, boundaryLoop[i], boundaryLoop[i + 1], boundaryLoop[i + 2]);
	}
}

void fetchCentroidSigmas(Mesh* mesh, vector<Eigen::Vector3i> triangles, vector<float>& centroidSigmas, unordered_map<int,float> sigmas)
{
	for (int i = 0; i < triangles.size(); i++)
	{
		int v1 = triangles[i](0);
		int v2 = triangles[i](1);
		int v3 = triangles[i](2);

		centroidSigmas.push_back((sigmas[v1] + sigmas[v2] + sigmas[v3]) / 3.0);
	}
}

void calculateSigmas(Mesh* mesh, const vector<int>& boundaryLoop, unordered_map<int, float>& sigmas)
{
	for (int i = 0; i < boundaryLoop.size(); i++)
	{
		sigmas[boundaryLoop[i]] = 0;

		for (int j = 0; j < mesh->verts[boundaryLoop[i]]->vertList.size(); j++)
		{
			float* adjc = mesh->verts[mesh->verts[boundaryLoop[i]]->vertList[j]]->coords;
			float* vc = mesh->verts[boundaryLoop[i]]->coords;

			sigmas[boundaryLoop[i]] += sqrt(pow(vc[0] - adjc[0], 2) + pow(vc[1] - adjc[1], 2) + pow(vc[2] - adjc[2], 2));
		}

		sigmas[boundaryLoop[i]] /= mesh->verts[boundaryLoop[i]]->vertList.size();
	}
}

void identifyDividedTriangles(Mesh* mesh, vector<Eigen::Vector3i> triangles,  vector<Eigen::Vector3d> centroids, unordered_map<int, float>sigmas, vector<float> centroidSigmas, vector<int>& trianglesToBeDivided)
{
	for (int i = 0; i < triangles.size(); i++)
	{
		float sigma_vc = centroidSigmas[i];
		Eigen::Vector3d vc = centroids[i];

		for (int m = 0; m < 3; m++)
		{
			Eigen::Vector3d vm;
			vm << mesh->verts[triangles[i][m]]->coords[0], mesh->verts[triangles[i][m]]->coords[1], mesh->verts[triangles[i][m]]->coords[2];
			float sigma_vm = sigmas[triangles[i][m]];

			if (ALPHA * (vc - vm).norm() > sigma_vc && ALPHA * (vc - vm).norm() > sigma_vm)
			{
				trianglesToBeDivided.push_back(i);
				break;
			}
		}
	}
}

void relaxEdge(Mesh* mesh, vector<Eigen::Vector3i>& triangles, int v1, int v2, vector<Eigen::Vector3d>& centroids, unordered_map<int, float>& sigmas, vector<float>& centroidSigmas)
{
	vector<int> trs;

	for (int i = 0; i < triangles.size(); i++)
	{
		unordered_set<int> tr;
		tr.insert(triangles[i](0));
		tr.insert(triangles[i](1));
		tr.insert(triangles[i](2));

		if (tr.count(v1) > 0 && tr.count(v2) > 0)
		{
			trs.push_back(i);
		}
	}

	if (trs.size() != 2)
		return;

	int op1 = -1;
	for (int i = 0; i < 3; i++)
	{
		if (triangles[trs[0]](i) != v1 && triangles[trs[0]](i) != v2)
		{
			op1 = triangles[trs[0]](i);
			break;
		}
	}

	int op2 = -1;
	for (int i = 0; i < 3; i++)
	{
		if (triangles[trs[1]](i) != v1 && triangles[trs[1]](i) != v2)
		{
			op2 = triangles[trs[1]](i);
			break;
		}
	}

	if (op1 == -1 || op2 == -1)
		return;

	if (circumSphereCheck(mesh, v1, v2, op1, op2) || circumSphereCheck(mesh, v1, v2, op2, op1))
	{
		if (swapCheck(triangles, v1, v2, op2, op1))
		{
			triangles[trs[0]] << op1, op2, v1;
			triangles[trs[1]] << op1, op2, v2;

			centroids[trs[0]] << (mesh->verts[op1]->coords[0] + mesh->verts[op2]->coords[0] + mesh->verts[v1]->coords[0]) / 3.0,
				(mesh->verts[op1]->coords[1] + mesh->verts[op2]->coords[1] + mesh->verts[v1]->coords[1]) / 3.0,
				(mesh->verts[op1]->coords[2] + mesh->verts[op2]->coords[2] + mesh->verts[v1]->coords[2]) / 3.0;

			centroids[trs[1]] << (mesh->verts[op1]->coords[0] + mesh->verts[op2]->coords[0] + mesh->verts[v2]->coords[0]) / 3.0,
				(mesh->verts[op1]->coords[1] + mesh->verts[op2]->coords[1] + mesh->verts[v2]->coords[1]) / 3.0,
				(mesh->verts[op1]->coords[2] + mesh->verts[op2]->coords[2] + mesh->verts[v2]->coords[2]) / 3.0;

			centroidSigmas[trs[0]] = (sigmas[op1] + sigmas[op2] + sigmas[v1]) / 3.0;
			centroidSigmas[trs[1]] = (sigmas[op1] + sigmas[op2] + sigmas[v2]) / 3.0;
		}
	}
}

void divideTriangles(Mesh* mesh, vector<Eigen::Vector3i>& triangles, const vector<int>& trianglesToBeDivided, vector<Eigen::Vector3d>& centroids, unordered_map<int, float>& sigmas, vector<float>& centroidSigmas, vector<pair<int, int>>& edgesToBeRelaxed)
{
	for (int i = 0; i < trianglesToBeDivided.size(); i++)
	{
		float* coords = (float*)malloc(sizeof(float) * 3);
		coords[0] = centroids[trianglesToBeDivided[i]](0);
		coords[1] = centroids[trianglesToBeDivided[i]](1);
		coords[2] = centroids[trianglesToBeDivided[i]](2);
		
		Vertex* centroid = new Vertex(mesh->verts.size(), coords);
		mesh->verts.push_back(centroid);

		edgesToBeRelaxed.push_back(pair<int, int>(triangles[trianglesToBeDivided[i]](0), triangles[trianglesToBeDivided[i]](1)));
		edgesToBeRelaxed.push_back(pair<int, int>(triangles[trianglesToBeDivided[i]](1), triangles[trianglesToBeDivided[i]](2)));
		edgesToBeRelaxed.push_back(pair<int, int>(triangles[trianglesToBeDivided[i]](2), triangles[trianglesToBeDivided[i]](0)));

		Eigen::Vector3i tr1;
		tr1 << mesh->verts.size() - 1, triangles[trianglesToBeDivided[i]](1), triangles[trianglesToBeDivided[i]](2);
		Eigen::Vector3i tr2;
		tr2 << triangles[trianglesToBeDivided[i]](0), mesh->verts.size() - 1, triangles[trianglesToBeDivided[i]](2);
		Eigen::Vector3i tr3;
		tr3 << triangles[trianglesToBeDivided[i]](0), triangles[trianglesToBeDivided[i]](1), mesh->verts.size() - 1;
		triangles[trianglesToBeDivided[i]] = tr1;
		triangles.push_back(tr2);
		triangles.push_back(tr3);

		sigmas[mesh->verts.size() - 1] = centroidSigmas[trianglesToBeDivided[i]];

		//calculate centroids & centroidsigmas
		Eigen::Vector3d c1;
		c1 << (mesh->verts[tr1(0)]->coords[0] + mesh->verts[tr1(1)]->coords[0] + mesh->verts[tr1(2)]->coords[0]) / 3.0,
			(mesh->verts[tr1(0)]->coords[1] + mesh->verts[tr1(1)]->coords[1] + mesh->verts[tr1(2)]->coords[1]) / 3.0,
			(mesh->verts[tr1(0)]->coords[2] + mesh->verts[tr1(1)]->coords[2] + mesh->verts[tr1(2)]->coords[2]) / 3.0;
		Eigen::Vector3d c2;
		c2 << (mesh->verts[tr2(0)]->coords[0] + mesh->verts[tr2(1)]->coords[0] + mesh->verts[tr2(2)]->coords[0]) / 3.0,
			(mesh->verts[tr2(0)]->coords[1] + mesh->verts[tr2(1)]->coords[1] + mesh->verts[tr2(2)]->coords[1]) / 3.0,
			(mesh->verts[tr2(0)]->coords[2] + mesh->verts[tr2(1)]->coords[2] + mesh->verts[tr2(2)]->coords[2]) / 3.0;
		Eigen::Vector3d c3;
		c3 << (mesh->verts[tr3(0)]->coords[0] + mesh->verts[tr3(1)]->coords[0] + mesh->verts[tr3(2)]->coords[0]) / 3.0,
			(mesh->verts[tr3(0)]->coords[1] + mesh->verts[tr3(1)]->coords[1] + mesh->verts[tr3(2)]->coords[1]) / 3.0,
			(mesh->verts[tr3(0)]->coords[2] + mesh->verts[tr3(1)]->coords[2] + mesh->verts[tr3(2)]->coords[2]) / 3.0;

		centroidSigmas[trianglesToBeDivided[i]] = (sigmas[tr1(0)] + sigmas[tr1(1)] + sigmas[tr1(2)]) / 3.0;
		centroidSigmas.push_back((sigmas[tr2(0)] + sigmas[tr2(1)] + sigmas[tr2(2)]) / 3.0);
		centroidSigmas.push_back((sigmas[tr3(0)] + sigmas[tr3(1)] + sigmas[tr3(2)]) / 3.0);

		centroids[trianglesToBeDivided[i]] = c1;
		centroids.push_back(c2);
		centroids.push_back(c3);

		
	}
}

bool circumSphereCheck(Mesh* mesh, int v1_, int v2_, int v3_, int op_)
{
	Eigen::Vector3d v1;
	v1 << mesh->verts[v1_]->coords[0], mesh->verts[v1_]->coords[1], mesh->verts[v1_]->coords[2];
	Eigen::Vector3d v2;
	v2 << mesh->verts[v2_]->coords[0], mesh->verts[v2_]->coords[1], mesh->verts[v2_]->coords[2];
	Eigen::Vector3d v3;
	v3 << mesh->verts[v3_]->coords[0], mesh->verts[v3_]->coords[1], mesh->verts[v3_]->coords[2];
	Eigen::Vector3d op;
	op << mesh->verts[op_]->coords[0], mesh->verts[op_]->coords[1], mesh->verts[op_]->coords[2];

	if ((v1 - v2).norm() < (v3 - op).norm())
		return false;

	Eigen::Matrix<double, 2, 3> projection; // 3D to 2D projection
	Eigen::Vector3d v10 = (v2 - v1).normalized();
	Eigen::Vector3d v = v3 - v1;
	Eigen::Vector3d n = v.cross(v10);
	Eigen::Vector3d v20 = v10.cross(n).normalized();

	projection.row(0) = v10;
	projection.row(1) = v20;

	Eigen::Vector2d A;
	A.fill(0);

	Eigen::Matrix<double, 2, 1> B = projection * (v2 - v1);
	Eigen::Matrix<double, 2, 1> C = projection * (v3 - v1);
	Eigen::Matrix<double, 2, 1> v11 = projection * (op - v1);

	Eigen::Matrix4d M;
	M(0, 0) = v11.squaredNorm();
	M(1, 0) = A.squaredNorm();
	M(2, 0) = B.col(0).squaredNorm();
	M(3, 0) = C.col(0).squaredNorm();

	M(0, 1) = v11(0);
	M(1, 1) = A(0);
	M(2, 1) = B(0);
	M(3, 1) = C(0);

	M(0, 2) = v11(1);
	M(1, 2) = A(1);
	M(2, 2) = B(1);
	M(3, 2) = C(1);

	M(0, 3) = 1;
	M(1, 3) = 1;
	M(2, 3) = 1;
	M(3, 3) = 1;

	if (M.determinant() > 0)
		return false;
	else
		return true;
}

bool swapCheck(const vector<Eigen::Vector3i>& triangles, int iv1, int iv2, int op2, int op1)
{
	unordered_set<int> tr1, tr2;
	tr1.insert(iv1);
	tr1.insert(op2);
	tr1.insert(op1);
	tr2.insert(iv2);
	tr2.insert(op1);
	tr2.insert(op2);

	for (int i = 0; i < triangles.size(); i++)
	{
		if (tr1.count(triangles[i](0)) > 0 && tr1.count(triangles[i](1)) > 0 && tr1.count(triangles[i](2)) > 0)
			return false;

		if (tr2.count(triangles[i](0)) > 0 && tr2.count(triangles[i](1)) > 0 && tr2.count(triangles[i](2)) > 0)
			return false;
	}

	return true;
}


void swap(Mesh* mesh, vector<Eigen::Vector3i>& triangles, int tr1, int op1, int tr2, int op2, int iv1, int iv2, vector<Eigen::Vector3d>& centroids, unordered_map<int, float>& sigmas, vector<float>& centroidSigmas)
{
	if (swapCheck(triangles, iv1, iv2, op2, op1))
	{
		triangles[tr1] << iv1, op2, op1;
		triangles[tr2] << iv2, op1, op2;

		centroids[tr1] << (mesh->verts[op1]->coords[0] + mesh->verts[op2]->coords[0] + mesh->verts[iv1]->coords[0]) / 3.0,
			(mesh->verts[op1]->coords[1] + mesh->verts[op2]->coords[1] + mesh->verts[iv1]->coords[1]) / 3.0,
			(mesh->verts[op1]->coords[2] + mesh->verts[op2]->coords[2] + mesh->verts[iv1]->coords[2]) / 3.0;

		centroids[tr2] << (mesh->verts[op1]->coords[0] + mesh->verts[op2]->coords[0] + mesh->verts[iv2]->coords[0]) / 3.0,
			(mesh->verts[op1]->coords[1] + mesh->verts[op2]->coords[1] + mesh->verts[iv2]->coords[1]) / 3.0,
			(mesh->verts[op1]->coords[2] + mesh->verts[op2]->coords[2] + mesh->verts[iv2]->coords[2]) / 3.0;

		centroidSigmas[tr1] = (sigmas[op1] + sigmas[op2] + sigmas[iv1]) / 3.0;
		centroidSigmas[tr2] = (sigmas[op1] + sigmas[op2] + sigmas[iv2]) / 3.0;
	}
}

struct swapS
{
	int tr1;
	int op1;
	int tr2;
	int op2;
	int iv1;
	int iv2;
};

bool edgeCheck(Mesh* mesh, vector<Eigen::Vector3i>& triangles, unordered_map<edge, int>& edgesToTriangles, int v1, int v2, int v3, int i, unordered_set<int>& swapSet, vector<swapS>& toBeSwapped)
{
	int a, b;

	if ((a = edgesToTriangles.count({ triangles[i](v1),  triangles[i](v2) })) == 0 && (b = edgesToTriangles.count({ triangles[i](v2),  triangles[i](v1) })) == 0)
	{
		edgesToTriangles[{ triangles[i](v1), triangles[i](v2) }] = i;
	}
	else
	{
		int neighbor, opposite  = -1;

		neighbor = (a > 0) ? edgesToTriangles[{ triangles[i](v1), triangles[i](v2) }] : edgesToTriangles[{ triangles[i](v2), triangles[i](v1) }];

		if (swapSet.count(neighbor))
			return false;

		for (int k = 0; k < 3; k++)
		{
			if (triangles[neighbor](k) != triangles[i](v1) && triangles[neighbor](k) != triangles[i](v2))
			{
				opposite = triangles[neighbor](k);
				break;
			}
		}
		
		if (opposite == -1)
		{
			return false;
		}
		
		if (circumSphereCheck(mesh, triangles[i](v1), triangles[i](v2), triangles[i](v3), opposite) || circumSphereCheck(mesh, triangles[i](v1), triangles[i](v2), opposite, triangles[i](v3)))
		{
			swapSet.insert(i);
			swapSet.insert(neighbor);
			toBeSwapped.push_back({ i, triangles[i](v3), neighbor, opposite, triangles[i](v1), triangles[i](v2) });
			return true;
		}
	}

	return false;
}

void relax(Mesh* mesh, vector<Eigen::Vector3i>& triangles, int& numOfSwap, int& sameCounter, vector<Eigen::Vector3d>& centroids, unordered_map<int, float>& sigmas, vector<float>& centroidSigmas)
{
	unordered_map<edge, int> edgesToTriangles;
	unordered_set<int> swapSet;
	vector<swapS> toBeSwapped;

	for (int i = 0; i < triangles.size(); i++)
	{
		if (swapSet.count(i))
			continue;

		edge e1{ triangles[i](0),  triangles[i](1) };
		edge e2{ triangles[i](1),  triangles[i](2) };
		edge e3{ triangles[i](2),  triangles[i](0) };
		
		if (edgeCheck(mesh, triangles, edgesToTriangles, 0, 1, 2, i, swapSet, toBeSwapped)) continue;
		if (edgeCheck(mesh, triangles, edgesToTriangles, 1, 2, 0, i, swapSet, toBeSwapped)) continue;
		if (edgeCheck(mesh, triangles, edgesToTriangles, 2, 0, 1, i, swapSet, toBeSwapped)) continue;
	}

	for (int i = 0; i < toBeSwapped.size(); i++)
	{
		swap(mesh, triangles, toBeSwapped[i].tr1, toBeSwapped[i].op1, toBeSwapped[i].tr2, toBeSwapped[i].op2, toBeSwapped[i].iv1, toBeSwapped[i].iv2, centroids, sigmas, centroidSigmas);
	}

	if (toBeSwapped.size() == numOfSwap)
		sameCounter++;

	numOfSwap = toBeSwapped.size();
}

void refinement(Mesh* mesh, const vector<int>& boundaryLoop, vector<Eigen::Vector3i>& triangles, const vector<vector<int>>& ls, enum RMETHOD rmethod)
{
	vector<Eigen::Vector3d> centroids;
	triangles = trianglesToBeInserted(mesh, boundaryLoop, ls, centroids);

	if (rmethod == NO_REFINEMENT)
		return;

	unordered_map<int, float> sigmas;
	calculateSigmas(mesh, boundaryLoop, sigmas);

	vector<float> centroidSigmas;
	fetchCentroidSigmas(mesh, triangles, centroidSigmas, sigmas);

	while (true)
	{
		vector<int> trianglesToBeDivided;
		identifyDividedTriangles(mesh, triangles, centroids, sigmas, centroidSigmas, trianglesToBeDivided);

		if (trianglesToBeDivided.empty()) break;

		vector<pair<int, int>> edgesToBeRelaxed;

		divideTriangles(mesh, triangles, trianglesToBeDivided, centroids, sigmas, centroidSigmas, edgesToBeRelaxed);

		for (int i = 0; i < edgesToBeRelaxed.size(); i++)
		{
			relaxEdge(mesh, triangles, edgesToBeRelaxed[i].first, edgesToBeRelaxed[i].second, centroids, sigmas, centroidSigmas);
		}

		int numOfSwap = INT_MAX;
		int sameCounter = 0;

		while (numOfSwap && sameCounter < INFINITE_PREVENTER)
		{
			relax(mesh, triangles, numOfSwap, sameCounter, centroids, sigmas, centroidSigmas);
		}
	}
}

double weight(Eigen::Vector3d v1, Eigen::Vector3d v2, enum FMETHOD fmethod)
{
	switch (fmethod)
	{
		case UNIFORM:
			return 1;
		case SCALE:
			return (v1 - v2).norm();
	}
	
}

double weightHarmonic(Mesh* mesh, int v1, int v2)
{
	vector<int> tris;

	for (int i = 0; i < mesh->tris.size(); i++)
	{
		unordered_set<int> verts;

		verts.insert(mesh->tris[i]->v1i);
		verts.insert(mesh->tris[i]->v2i);
		verts.insert(mesh->tris[i]->v3i);

		if (verts.count(v1) > 0 && verts.count(v2) > 0)
			tris.push_back(i);
	}

	if (tris.size() != 2)
		return 1;

	Eigen::Vector3i tri1, tri2;
	tri1 << mesh->tris[tris[0]]->v1i, mesh->tris[tris[0]]->v2i, mesh->tris[tris[0]]->v3i;
	tri2 << mesh->tris[tris[1]]->v1i, mesh->tris[tris[1]]->v2i, mesh->tris[tris[1]]->v3i;

	int op1, op2;

	for (int i = 0; i < 3; i++)
	{
		if (tri1(i) != v1 && tri1(i) != v2)
			op1 = tri1(i);
		if (tri2(i) != v1 && tri2(i) != v2)
			op2 = tri2(i);
	}

	Eigen::Vector3d v11, v22, op11, op22;
	v11 << mesh->verts[v1]->coords[0], mesh->verts[v1]->coords[1], mesh->verts[v1]->coords[2];
	v22 << mesh->verts[v2]->coords[0], mesh->verts[v2]->coords[1], mesh->verts[v2]->coords[2];
	op11 << mesh->verts[op1]->coords[0], mesh->verts[op1]->coords[1], mesh->verts[op1]->coords[2];
	op22 << mesh->verts[op2]->coords[0], mesh->verts[op2]->coords[1], mesh->verts[op2]->coords[2];

	double w1 = 1/tan(acos((v11 - op11).dot((v22 - op11)) / ((v11 - op11).norm() * (v22 - op11).norm())));
	double w2 = 1 / tan(acos((v11 - op22).dot((v22 - op22)) / ((v11 - op22).norm() * (v22 - op22).norm())));

	return w1 + w2;
}

void fairUniformWhole(Mesh* mesh, enum FMETHOD fmethod)
{
	for (int v = 0; v < mesh->verts.size(); v++)
	{
		Eigen::Vector3d vv;
		vv << mesh->verts[v]->coords[0], mesh->verts[v]->coords[1], mesh->verts[v]->coords[2];

		Eigen::Vector3d Uw = -vv;

		double wv = 0;
		Eigen::Vector3d sum_wv(0, 0, 0);

		for (int k = 0; k < mesh->verts[v]->vertList.size(); k++)
		{
			Eigen::Vector3d vi;
			vi << mesh->verts[mesh->verts[v]->vertList[k]]->coords[0], mesh->verts[mesh->verts[v]->vertList[k]]->coords[1], mesh->verts[mesh->verts[v]->vertList[k]]->coords[2];

			double w = weight(vv, vi, (fmethod == UNIFORM_WHOLE) ? UNIFORM : SCALE);

			wv += w;
			sum_wv += vi * w;
		}

		Uw += sum_wv / wv;

		vv += Uw;
		mesh->verts[v]->coords[0] = vv(0);
		mesh->verts[v]->coords[1] = vv(1);
		mesh->verts[v]->coords[2] = vv(2);
	}

}


void fairing(Mesh* mesh, const vector<Eigen::Vector3i>& triangles, enum FMETHOD fmethod)
{
	switch (fmethod)
	{
	case NO_FAIR:
		return;
	case UNIFORM_WHOLE:
		fairUniformWhole(mesh, UNIFORM_WHOLE);
		return;
	case SCALE_WHOLE:
		fairUniformWhole(mesh, SCALE_WHOLE);
		return;
	default:
		break;
	}

	unordered_set<int> verticesProcessed;

	for (int i = 0; i < triangles.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int v = triangles[i](j);

			if (verticesProcessed.count(v) > 0)
				continue;

			Eigen::Vector3d vv;
			vv << mesh->verts[v]->coords[0], mesh->verts[v]->coords[1], mesh->verts[v]->coords[2];

			Eigen::Vector3d Uw = -vv;

			double wv = 0;
			Eigen::Vector3d sum_wv(0,0,0);

			for (int k = 0; k < mesh->verts[v]->vertList.size(); k++)
			{
				Eigen::Vector3d vi;
				vi << mesh->verts[mesh->verts[v]->vertList[k]]->coords[0], mesh->verts[mesh->verts[v]->vertList[k]]->coords[1], mesh->verts[mesh->verts[v]->vertList[k]]->coords[2];

				double w;

				if (fmethod != HARMONIC)
				{
					w = weight(vv, vi, fmethod);
				}
				else
				{
					w = weightHarmonic(mesh, v, mesh->verts[v]->vertList[k]);
				}

				wv += w;
				sum_wv += vi * w;
			}

			Uw += sum_wv / wv;

			vv += Uw;
			mesh->verts[v]->coords[0] = vv(0);
			mesh->verts[v]->coords[1] = vv(1);
			mesh->verts[v]->coords[2] = vv(2);
			
			verticesProcessed.insert(v);
		}
	}

}

void holyFillerHelper(Mesh* mesh, const vector<int>& boundaryLoop, vector<Eigen::Vector3i>& filled, enum METHOD method, enum RMETHOD rmethod, enum FMETHOD fmethod, int numFairIt)
{
	vector<Eigen::Vector3i> triangles;
	vector<vector<double>> A;
	vector<vector<int>> ls;
	
	for (int i = boundaryLoop.size() - 1; i > 0; i--)
	{	
		vector<double> vec(i);
		A.push_back(vec);
		
		if (i >= boundaryLoop.size() - 2)
		{
			vector<int> v;
			ls.push_back(v);
		}
		else
		{
			vector<int> v(i);
			ls.push_back(v);
		}
	}

	if (method == ANGLE)
	{
		vector<vector<Eigen::Vector3d>> edgeTriangleNormals;
		vector<vector<double>> dotPs;

		findEdgeTriangleNormals(mesh, edgeTriangleNormals, boundaryLoop);

		for (int i = boundaryLoop.size() - 1; i > 0; i--)
		{
			vector<double> vec(i, 1);
			dotPs.push_back(vec);
		}

		for (int i = 3; i < boundaryLoop.size(); i++)
		{
			for (int j = 0; j < boundaryLoop.size() - i; j++)
			{
				Eigen::Vector3d optNormal;
				int KoPT = INT_MIN;
				double minArea = DBL_MAX;
				double maxDOTP = DBL_MIN;

				for (int k = 0; k < i - 1; k++)
				{
					double area = trArea(mesh->verts, boundaryLoop[j], boundaryLoop[j + k + 1], boundaryLoop[j + i]) + A[i - k - 2][j + k + 1] + A[k][j];
					Eigen::Vector3d normal = trNormal(mesh->verts, boundaryLoop[j], boundaryLoop[j + k + 1], boundaryLoop[j + i]);
					double dotPro = min(normal.dot(edgeTriangleNormals[k][j]), normal.dot(edgeTriangleNormals[i - k - 2][j + k + 1]));

					if (j == 0 && i == boundaryLoop.size() - 1)
						dotPro = min(dotPro, normal.dot(edgeTriangleNormals[0][boundaryLoop.size() - 1]));

					dotPro = min(min(dotPro, dotPs[k][j]), dotPs[i - k - 2][j + k + 1]);

					if ((maxDOTP == dotPro && area < minArea) || maxDOTP < dotPro)
					{
						optNormal = normal;
						maxDOTP = dotPro;
						KoPT = k;
						minArea = area;
					}
				}

				edgeTriangleNormals[i - 1][j] = optNormal;
				ls[i - 1][j] = KoPT + j + 1;
				A[i - 1][j] = minArea;
				dotPs[i - 1][j] = maxDOTP;
			}
		}
	}
	else //// AREA
	{
		for (int i = 3; i < boundaryLoop.size(); i++)
		{
			for (int j = 0; j < boundaryLoop.size() - i; j++)
			{
				double minArea = DBL_MAX;
				int KoPT = INT_MIN;
				
				for (int k = 0; k < i - 1; k++)
				{
					double area = trArea(mesh->verts, boundaryLoop[j], boundaryLoop[k + j + 1], boundaryLoop[j + i]) + A[k][j] + A[i - k - 2][k + j + 1];
					
					if (area < minArea)
					{
						KoPT = k;
						minArea = area;
					}
				}
				
				ls[i - 1][j] = KoPT + j + 1;
				A[i - 1][j] = minArea;
			}
		}
	}

	refinement(mesh, boundaryLoop, triangles, ls, rmethod);
	
	insertTriangles(mesh, triangles);
	
	for(int k = 0; k < numFairIt; k++) 
		fairing(mesh, triangles, fmethod);

	filled.insert(filled.end(), triangles.begin(), triangles.end());
}

vector<Eigen::Vector3i> holyFiller(Mesh* mesh, enum METHOD method, enum RMETHOD rmethod, enum FMETHOD fmethod, int numFairIt)
{
	vector<vector<int>> boundaryLoops;
	boundaryLoopDetector(mesh->tris, boundaryLoops);
	vector<Eigen::Vector3i> filled;
	for (auto boundaryLoop : boundaryLoops)
	{
		holyFillerHelper(mesh, boundaryLoop, filled, method, rmethod, fmethod, numFairIt);
	}

	return filled;
}
