#include "Mesh.h"

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

void insertTriangles(Mesh* mesh, const vector<int>& boundaryLoop, const vector<vector<int>>& lambdas, vector<int>& filled)
{
	vector<pair<int, int>> sections;

	sections.push_back(pair<int, int>{ 0, boundaryLoop.size() - 1 });

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

		filled.push_back(boundaryLoop[section.first]);
		filled.push_back(boundaryLoop[k]);
		filled.push_back(boundaryLoop[section.second]);
		mesh->tris.push_back(new Triangle(mesh->tris.size(), boundaryLoop[section.first], boundaryLoop[k], boundaryLoop[section.second]));

		if (k - section.first > 1)
		{
			sections.push_back(pair<int, int>{ section.first, k });
		}

		if (section.second - k > 1)
		{
			sections.push_back(pair<int, int>{ k, section.second });
		}
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

void holyFillerHelper(Mesh* mesh, const vector<int>& boundaryLoop, vector<int>& filled, enum METHOD method)
{
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

	insertTriangles(mesh, boundaryLoop, ls, filled);
}


vector<int> holyFiller(Mesh* mesh, enum METHOD method)
{
	vector<int> filled;
	vector<vector<int>> boundaryLoops;
	boundaryLoopDetector(mesh->tris, boundaryLoops);

	for (auto boundaryLoop : boundaryLoops)
	{
		holyFillerHelper(mesh, boundaryLoop, filled, method);
	}

	return filled;
}
