#include "Utils.hpp"
#include "PolyhedronMesh.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <queue>
#include <algorithm>
#include <map>
#include <set>
#include <Eigen/Dense>
#include "UCDUtilities.hpp"


using namespace std;
using namespace PolyhedronMesh;

void normalize(vertex& v) { //O(1)
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len < 1e-8) return; 
    v.x /= len;
    v.y /= len;
    v.z /= len;
}

void buildTetrahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) { // O(v+e+f), O(1) con dati piccoli
	vertices.clear();
	edges.clear();
	faces.clear();
	
	vertices = {
		{0, 1, 1, 1},
		{1, -1, -1, 1},
		{2, -1, 1, -1},
		{3, 1, -1, -1}
	};
	
	edges = {
		{0, 0, 1}, {1, 0, 2}, {2, 0, 3},
		{3, 1, 2}, {4, 1, 3}, {5, 2, 3}
	};
	
	faces = {
		{0, {0, 1, 2}, {0, 3, 1}},  
		{1, {0, 2, 3}, {1, 5, 2}},   
		{2, {0, 3, 1}, {2, 4, 0}},   
		{3, {1, 2, 3}, {3, 5, 4}}    
	};
	
	polyhedron.id = 0;
	polyhedron.vertex_ids.clear();
	for (size_t i = 0; i < vertices.size(); ++i) {
		polyhedron.vertex_ids.push_back(i);
	}

	polyhedron.edge_ids.clear();
	for (size_t i = 0; i < edges.size(); ++i) {
		polyhedron.edge_ids.push_back(i);
	}

	polyhedron.face_ids.clear();
	for (size_t i = 0; i < faces.size(); ++i) {
		polyhedron.face_ids.push_back(i);
	}
}

void buildEsahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) { // O(v+e+f), O(1) con dati piccoli
	vertices.clear();
	edges.clear();
	faces.clear();
	
	vertices = {
		{0, -1, -1, -1},
		{1,  1, -1, -1},
		{2,  1,  1, -1},
		{3, -1,  1, -1},
		{4, -1, -1,  1},
		{5,  1, -1,  1},
		{6,  1,  1,  1},
		{7, -1,  1,  1}
	};
	
	edges = {
		{0, 0, 1}, {1, 1, 2}, {2, 2, 3}, {3, 3, 0},
        {4, 4, 5}, {5, 5, 6}, {6, 6, 7}, {7, 7, 4}, 
        {8, 0, 4}, {9, 1, 5}, {10, 2, 6}, {11, 3, 7}
	};
		
	faces = {
		{0, {0, 1, 2, 3},  {0, 1, 2, 3}},   
		{1, {4, 5, 6, 7},  {4, 5, 6, 7}},   
		{2, {0, 1, 5, 4},  {0, 9, 4, 8}},   
		{3, {1, 2, 6, 5},  {1,10, 5, 9}},   
		{4, {2, 3, 7, 6},  {2,11, 6,10}},   
		{5, {3, 0, 4, 7},  {3, 8, 7,11}}    
	};

	
	polyhedron.id = 1;
    polyhedron.vertex_ids.clear();
	for (size_t i = 0; i < vertices.size(); ++i) {
		polyhedron.vertex_ids.push_back(i);
	}

	polyhedron.edge_ids.clear();
	for (size_t i = 0; i < edges.size(); ++i) {
		polyhedron.edge_ids.push_back(i);
	}

	polyhedron.face_ids.clear();
	for (size_t i = 0; i < faces.size(); ++i) {
		polyhedron.face_ids.push_back(i);
	}
}

void buildOctahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
	
	vertices.clear();
	edges.clear();
	faces.clear();
	
	vertices = {
        {0,  1,  0,  0},
        {1, -1,  0,  0},  
        {2,  0,  1,  0},  
        {3,  0, -1,  0},  
        {4,  0,  0,  1}, 
        {5,  0,  0, -1}   
    };
	
	edges = {
        {0, 0, 2}, {1, 2, 1}, {2, 1, 3}, {3, 3, 0}, 
        {4, 0, 4}, {5, 2, 4}, {6, 1, 4}, {7, 3, 4}, 
        {8, 0, 5}, {9, 2, 5}, {10, 1, 5}, {11, 3, 5} 
    };
	
	faces = {
		{0, {0, 2, 4}, {0, 5, 4}},     
		{1, {2, 1, 4}, {1, 6, 5}},     
		{2, {1, 3, 4}, {2, 7, 6}},     
		{3, {3, 0, 4}, {3, 4, 7}},     
		{4, {0, 5, 2}, {8, 9, 0}},     
		{5, {2, 5, 1}, {9,10, 1}},     
		{6, {1, 5, 3}, {10,11, 2}},    
		{7, {3, 5, 0}, {11, 8, 3}}     
	};

	polyhedron.id = 2;
    polyhedron.vertex_ids.clear();
	for (size_t i = 0; i < vertices.size(); ++i) {
		polyhedron.vertex_ids.push_back(i);
	}

	polyhedron.edge_ids.clear();
	for (size_t i = 0; i < edges.size(); ++i) {
		polyhedron.edge_ids.push_back(i);
	}

	polyhedron.face_ids.clear();
	for (size_t i = 0; i < faces.size(); ++i) {
		polyhedron.face_ids.push_back(i);
	}
	
}

void buildDodecahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
    vertices.clear();
    edges.clear();
    faces.clear();

    const double phi = (1.0 + sqrt(5.0)) / 2.0;
    const double a = 1.0;
    const double b = 1.0 / phi;
    const double c = phi;

    vector<Eigen::Vector3d> raw = {
        { a,  a,  a}, { a,  a, -a}, { a, -a,  a}, { a, -a, -a},
        {-a,  a,  a}, {-a,  a, -a}, {-a, -a,  a}, {-a, -a, -a},
        { 0,  b,  c}, { 0,  b, -c}, { 0, -b,  c}, { 0, -b, -c},
        { b,  c, 0}, {-b,  c,  0}, { b, -c,  0}, {-b, -c, 0},
        { c,  0,  b}, {-c,  0,  b}, { c,  0, -b}, {-c,  0, -b}
    };

    for (size_t i = 0; i < raw.size(); ++i) {
        raw[i].normalize();  
        vertices.push_back({(int)i, raw[i][0], raw[i][1], raw[i][2]});
    }

    vector<vector<int>> faceVertices = {
        {0, 8, 4, 13, 12},
        {0, 12, 1, 18, 16},
        {0, 16, 2, 10, 8},
        {8, 10, 6, 17, 4},
        {12, 13, 5, 9, 1},
        {16, 18, 3, 14, 2},
        {4, 17, 19, 5, 13},
        {1, 9, 11, 3, 18},
        {2, 14, 15, 6, 10},
        {3, 11, 7, 15, 14},
        {5, 19, 7, 11, 9},
        {6, 15, 7, 19, 17}
    };

    // mappa per garantire spigoli unici
    map<pair<int, int>, int> edgeMap;

    for (size_t i = 0; i < faceVertices.size(); ++i) {
        face f;
        f.id = static_cast<int>(i);
        f.vertex_ids = faceVertices[i];
        int n = f.vertex_ids.size();
        for (int j = 0; j < n; ++j) {
            int a = f.vertex_ids[j];
            int b = f.vertex_ids[(j + 1) % n];
            auto key = minmax(a, b);
            if (!edgeMap.count(key)) {
                int eid = static_cast<int>(edges.size());
                edgeMap[key] = eid;
                edges.push_back({eid, key.first, key.second});
            }
            f.edge_ids.push_back(edgeMap[key]);
        }
        faces.push_back(f);
    }

    polyhedron.id = 3;
    polyhedron.vertex_ids.clear();
    for (size_t i = 0; i < vertices.size(); ++i)
        polyhedron.vertex_ids.push_back(i);
    polyhedron.edge_ids.clear();
    for (size_t i = 0; i < edges.size(); ++i)
        polyhedron.edge_ids.push_back(i);
    polyhedron.face_ids.clear();
    for (size_t i = 0; i < faces.size(); ++i)
        polyhedron.face_ids.push_back(i);
}


void buildIcosahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
    vertices.clear();
	edges.clear();
	faces.clear();

    const double phi = (1.0 + sqrt(5.0)) / 2.0;

    vector<Eigen::Vector3d> points = {
        { 0,  1,  phi}, { 0, -1,  phi}, { phi, 0, 1}, { 1,  phi, 0}, {-1,  phi, 0}, {-phi, 0, 1},
        { 1, -phi, 0}, { phi, 0, -1}, { 0, 1, -phi}, {-phi, 0, -1}, {-1, -phi, 0}, { 0, -1, -phi}
    };

    for (size_t i = 0; i < points.size(); ++i) {
        Eigen::Vector3d p = points[i].normalized();
        vertices.push_back({static_cast<int>(i), p[0], p[1], p[2]});
    }

    vector<vector<int>> faceVertices = {
        {0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 5}, {0, 5, 1},
        {1, 6, 2}, {2, 6, 7}, {2, 7, 3}, {3, 7, 8}, {3, 8, 4},
        {4, 8, 9}, {4, 9, 5}, {5, 9,10}, {5,10, 1}, {1,10, 6},
        {11, 6,10}, {11, 7, 6}, {11, 8, 7}, {11, 9, 8}, {11,10,9}
    };

    for (size_t i = 0; i < faceVertices.size(); ++i) {
        face f;
        f.id = static_cast<int>(i);
        f.vertex_ids = faceVertices[i];
        faces.push_back(f);
    }

    set<pair<int, int>> edgeSet;
    for (auto& f : faces) {
        f.edge_ids.clear();
        for (int i = 0; i < 3; ++i) {
            int a = f.vertex_ids[i];
            int b = f.vertex_ids[(i + 1) % 3];
            auto key = std::minmax(a, b);
            if (!edgeSet.count(key)) {
                edgeSet.insert(key);
                edges.push_back({static_cast<int>(edges.size()), key.first, key.second});
            }
            for (const auto& e : edges) {
                if (minmax(e.origin, e.end) == key) {
                    f.edge_ids.push_back(e.id);
                    break;
                }
            }
        }
    }

    polyhedron.id = 4;
    polyhedron.vertex_ids.clear();
    for (size_t i = 0; i < vertices.size(); ++i)
        polyhedron.vertex_ids.push_back(i);
    polyhedron.edge_ids.clear();
    for (size_t i = 0; i < edges.size(); ++i)
        polyhedron.edge_ids.push_back(i);
    polyhedron.face_ids.clear();
    for (size_t i = 0; i < faces.size(); ++i)
        polyhedron.face_ids.push_back(i);
}

void buildPolyhedron(int p, int q, int b, int c, vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
    if (p < 3 || q < 3) {
        cout << "Valori non validi" << endl;
        return;
    }
	
	bool ClassI = (b == 0 && c > 0) || (b > 0 && c == 0);
    bool ClassII = (b == c && b > 0);
	
    if (!ClassI && !ClassII && (b != 0 || c != 0)) {
        cout << "Tipo di classe non supportato" << endl;
        return;
    }
    if (ClassI) {
        buildClassIGeodesic(p, q, b + c, vertices, edges, faces, polyhedron);
        return;
    }

    if (ClassII) {
        buildClassIIGeodesic(p, q, b, c, vertices, edges, faces, polyhedron);
        return;
    }

    if (p == 3 && q == 3 && b == 0 && c == 0) {
        buildTetrahedron(vertices, edges, faces, polyhedron);
        return;
    }

    if (p == 4 && q == 3 && b == 0 && c == 0) {
        buildEsahedron(vertices, edges, faces, polyhedron);
        return;
    }

    if (p == 3 && q == 4 && b == 0 && c == 0) {
        buildOctahedron(vertices, edges, faces, polyhedron);
        return;
    }

    if (p == 5 && q == 3 && b == 0 && c == 0) {
        buildDodecahedron(vertices, edges, faces, polyhedron);
        return;
    }

    if (p == 3 && q == 5 && b == 0 && c == 0) {
        buildIcosahedron(vertices, edges, faces, polyhedron);
        return;
    }
	
    cout << "Configurazione non supportata\n";
}


double distance(const vertex& a, const vertex& b) { 
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

bool sameVertex(const vertex& a, const vertex& b, double tolerance) { 
    return distance(a, b) <= tolerance;
}

int getOrAddVertex(double x, double y, double z, vector<vertex>& vertices, double tolerance) { 
    tolerance = 1e-3;
	for (const auto& existing : vertices) {
        if (sameVertex(existing, vertex{-1, x, y, z}, tolerance)) {
            return existing.id;
        }
    }
	
    int newId = vertices.size();
    vertex v{newId, x, y, z};
    vertices.push_back(v);
    return newId;
}

int getOrAddVertex(double x, double y, double z, vector<vertex>& vertices) { 
    return getOrAddVertex(x, y, z, vertices, 1e-4);
}

void projectVerticesOnUnitSphere(vector<vertex>&vertices) {
    for(auto& v : vertices){
        normalize(v);
    }
}

void buildClassIGeodesic(int p, int q, int b, vector<vertex>& vertices, vector<edge>& edges, vector<face>& faces, polyhedron& poly) { 
    vertices.clear();
    edges.clear();
    faces.clear();

    vector<vertex> baseVertices;
    vector<edge> baseEdges;
    vector<face> baseFaces;
    polyhedron basePoly;

    buildPolyhedron(p, q, 0, 0, baseVertices, baseEdges, baseFaces, basePoly);

    const int N = b;
    const double tolerance = 1e-6;

    for (const auto& f : baseFaces) {
        vertex A = baseVertices[f.vertex_ids[0]];
        vertex B = baseVertices[f.vertex_ids[1]];
        vertex C = baseVertices[f.vertex_ids[2]];

        Eigen::MatrixXi grid(N + 1, N + 1);
        grid.setConstant(-1);

        // crea vertici baricentrici
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N - i; ++j) {
                int k = N - i - j;

                double alpha = static_cast<double>(i) / N;
                double beta  = static_cast<double>(j) / N;
                double gamma = static_cast<double>(k) / N;

                double x = alpha * A.x + beta * B.x + gamma * C.x;
                double y = alpha * A.y + beta * B.y + gamma * C.y;
                double z = alpha * A.z + beta * B.z + gamma * C.z;

                vertex v = {-1, x, y, z};
                normalize(v);
                grid(i,j) = getOrAddVertex(v.x, v.y, v.z, vertices, tolerance);
            }
        }

        // triangola la griglia baricentrica
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N - i; ++j) {
                int v0 = grid(i,j);
                int v1 = grid(i+1,j);
                int v2 = grid(i,j+1);
                if (v0 >= 0 && v1 >= 0 && v2 >= 0)
                    faces.push_back({(int)faces.size(), {v0, v1, v2}, {}});

                if (i + j < N - 1) {
                    int v3 = grid(i+1,j+1);
                    if (v1 >= 0 && v2 >= 0 && v3 >= 0)
                        faces.push_back({(int)faces.size(), {v1, v3, v2}, {}});
                }
            }
        }
    }

    // crea spigoli unici e associa agli ID facce
    map<pair<int, int>, int> edgeMap;
    for (auto& f : faces) {
        for (int i = 0; i < 3; ++i) {
            int a = f.vertex_ids[i];
            int b = f.vertex_ids[(i + 1) % 3];
            auto key = minmax(a, b);
            if (!edgeMap.count(key)) {
                int id = edges.size();
                edges.push_back({id, key.first, key.second});
                edgeMap[key] = id;
            }
            f.edge_ids.push_back(edgeMap[key]);
        }
    }

	polyhedron newPoly;
	newPoly.id = basePoly.id + 1000;

	for (int i = 0; i < (int)vertices.size(); ++i)
		newPoly.vertex_ids.push_back(i);
	for (int i = 0; i < (int)edges.size(); ++i)
		newPoly.edge_ids.push_back(i);
	for (int i = 0; i < (int)faces.size(); ++i)
		newPoly.face_ids.push_back(i);

	poly = newPoly;

}

void buildClassIIGeodesic(int p, int q, int b, int c, vector<vertex>& vertices, vector<edge>& edges, vector<face>& faces, polyhedron& poly) {
    int N = b + c;
    const double tol = 1e-4;

    vector<vertex> baseVertices;
    vector<edge> baseEdges;
    vector<face> baseFaces;
    polyhedron base;

    buildPolyhedron(p, q, 0, 0, baseVertices, baseEdges, baseFaces, base);

    vertices.clear();
    edges.clear();
    faces.clear();

    map<tuple<pair<int,int>, pair<int,int>, pair<int,int>>, int> vertexMap;
    int nextFaceId = 0;

    for (const auto& f : baseFaces) {
        if (f.vertex_ids.size() != 3) continue;

        vector<int> ids = f.vertex_ids;
		
		Eigen::Vector3d center(0, 0, 0);
		for (int vid : ids) {
			center += Eigen::Vector3d(baseVertices[vid].x, baseVertices[vid].y, baseVertices[vid].z);
		}
		center.normalize();

		Eigen::Vector3d vecs[3];
		for (int i = 0; i < 3; ++i) {
			vecs[i] = Eigen::Vector3d(baseVertices[ids[i]].x, baseVertices[ids[i]].y, baseVertices[ids[i]].z) - center;
		}

		// calcola angoli rispetto al primo vettore (reference)
		double angles[3];
		for (int i = 0; i < 3; ++i) {
			Eigen::Vector3d ref = vecs[0];
			Eigen::Vector3d curr = vecs[i];
			Eigen::Vector3d cross = ref.cross(curr);
			double dot = ref.dot(curr);
			double sign = (cross.dot(center) > 0) ? 1.0 : -1.0;
			angles[i] = sign * acos(min(1.0, max(-1.0, dot / (ref.norm() * curr.norm()))));
		}

		// ordina ids secondo gli angoli
		for (int i = 0; i < 2; ++i) {
			for (int j = i + 1; j < 3; ++j) {
				if (angles[i] > angles[j]) {
					swap(angles[i], angles[j]);
					swap(ids[i], ids[j]);
				}
			}
		}
		
		int a = ids[0];
		int b = ids[1];
		int c = ids[2];

        const vertex& A = baseVertices[a];
        const vertex& B = baseVertices[b];
        const vertex& C = baseVertices[c];

        Eigen::Vector3d va(A.x, A.y, A.z);
        Eigen::Vector3d vb(B.x, B.y, B.z);
        Eigen::Vector3d vc(C.x, C.y, C.z);

        Eigen::MatrixXi grid = Eigen::MatrixXi::Constant(N + 1, N + 1, -1);

        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N - i; ++j) {
                int k = N - i - j;

                //rotazione della base baricentrica
				Eigen::Vector3d P = ((double)k / N) * va + ((double)i / N) * vb + ((double)j / N) * vc;
                P.normalize();

                vector<pair<int, int>> parts = {
                    {a, i}, {b, j}, {c, k}
                };
                sort(parts.begin(), parts.end());
                auto key = make_tuple(parts[0], parts[1], parts[2]);

                int vid;
                if (vertexMap.count(key)) {
                    vid = vertexMap[key];
                } else {
                    vid = getOrAddVertex(P[0], P[1], P[2], vertices, tol);
                    vertexMap[key] = vid;
                }

                grid(i,j) = vid;
            }
        }

        // crea triangoli
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N - i; ++j) {
                int v0 = grid(i,j);
                int v1 = grid(i+1,j);
                int v2 = grid(i,j+1);
                faces.push_back({nextFaceId++, {v0, v1, v2}, {}});

                if (i + j < N - 1) {
                    int v3 = grid(i+1,j+1);
                    faces.push_back({nextFaceId++, {v1, v3, v2}, {}});
                }
            }
        }
    }

    // costruzione spigoli unici
    map<pair<int, int>, int> edgeMap; //associa a una coppia ordinata di vertici(lo spigolo) una coppia ordinata di vertici
    for (auto& f : faces) {
        f.edge_ids.clear();
        for (int i = 0; i < 3; ++i) {
            int u = f.vertex_ids[i];
            int v = f.vertex_ids[(i + 1) % 3];
            auto key = minmax(u, v);
            if (!edgeMap.count(key)) {
                int eid = static_cast<int>(edges.size());
                edges.push_back({eid, key.first, key.second});
                edgeMap[key] = eid;
            }
            f.edge_ids.push_back(edgeMap[key]);
        }
    }

    poly.id = 1000 + N;
    poly.vertex_ids.clear();
    for (auto& v : vertices) poly.vertex_ids.push_back(v.id);
    poly.edge_ids.clear();
    for (auto& e : edges) poly.edge_ids.push_back(e.id);
    poly.face_ids.clear();
    for (auto& f : faces) poly.face_ids.push_back(f.id);
}

bool isFaceConsistent(const face& f, const vector<edge>& edges, double tolerance) {
    int n = f.vertex_ids.size();
    if (n < 3 || f.edge_ids.size() != n) return false;

    for (int i = 0; i < n; ++i) {
        int v_start = f.vertex_ids[i];
        int v_end   = f.vertex_ids[(i + 1) % n];
        int eid     = f.edge_ids[i];

        if (eid >= edges.size()) return false;

        const edge& e = edges[eid];

        // accetta anche direzione inversa
        bool ok = (e.origin == v_start && e.end == v_end) ||
                  (e.origin == v_end && e.end == v_start);

        if (!ok) {
            cerr << "Incoerenza nella faccia " << f.id << ": edge[" << eid << "] "
                      << "non collega " << v_start << " <-> " << v_end
                      << " (è " << e.origin << " → " << e.end << ")\n";
            return false;
        }
    }

    return true;
}

void findShortestPath(vector<vertex>& vertices, vector<edge>& edges, int startId, int endId) {
    int N = vertices.size();
    vector<double> dist(N, std::numeric_limits<double>::infinity());
    vector<int> prev(N, -1);
    vector<bool> visited(N, false);

    dist[startId] = 0.0;

    using P = pair<double, int>;
    priority_queue<P, vector<P>, greater<P>> pq;
    pq.push({0.0, startId});

    //costruzione grafo come lista adiacenza
    vector<vector<pair<int, double>>> adj(N);
    for (const auto& e : edges) {
        double len = e.length > 0 ? e.length : distance(vertices[e.origin], vertices[e.end]);
        adj[e.origin].emplace_back(e.end, len);
        adj[e.end].emplace_back(e.origin, len); // non orientato
    }

    //algoritmo di Dijkstra
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (visited[u]) continue;
        visited[u] = true;

        for (auto [v, w] : adj[u]) {
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.push({dist[v], v});
            }
        }
    }

    if (prev[endId] == -1) {
        cout << "Nessun cammino tra i vertici " << startId << " e " << endId << "\n";
        return;
    }

    // ricostruzione cammino
    vector<int> path;
    for (int at = endId; at != -1; at = prev[at])
        path.push_back(at);
    reverse(path.begin(), path.end());

    // resetta ShortPath
    for (auto& v : vertices) v.ShortPath = 0;
    for (auto& e : edges) e.ShortPath = 0;

    // marca i vertici
    for (int v : path)
        vertices[v].ShortPath = 1;

    // marca i lati e somma lunghezza
    int edgeCount = 0;
    double totalLength = 0.0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        int u = path[i], v = path[i + 1];

        for (auto& e : edges) {
            if ((e.origin == u && e.end == v) || (e.origin == v && e.end == u)) {
                e.ShortPath = 1;
                e.length = distance(vertices[e.origin], vertices[e.end]);
                totalLength += e.length;
                edgeCount++;
                break;
            }
        }
    }
    cout << "Numero di lati nel cammino minimo: " << edgeCount << endl;
    cout << "Somma delle lunghezze: " << totalLength << endl;
}

vector<vertex> calculateCentroids(const vector<vertex>& vertices, const vector<face>& faces) {
    vector<vertex> centroids;
    centroids.reserve(faces.size());
    for (const auto& f : faces) {
        double cx = 0, cy = 0, cz = 0;
        int n = f.vertex_ids.size();
        for (int vid : f.vertex_ids) {
            cx += vertices[vid].x;
            cy += vertices[vid].y;
            cz += vertices[vid].z;
        }
        cx /= n;
        cy /= n;
        cz /= n;

        vertex center(centroids.size(), cx, cy, cz);
        centroids.push_back(center);
    }
    return centroids;
}

void buildDualPolyhedron(const vector<vertex>& vertices,  const vector<edge>& edges, const vector<face>& faces, const polyhedron& original, vector<vertex>& dualVertices, vector<edge>& dualEdges, vector<face>& dualFaces,polyhedron& dualPoly) {
    dualVertices.clear();
    dualEdges.clear();
    dualFaces.clear();

    dualVertices = calculateCentroids(vertices, faces);
	for (auto& v : dualVertices) {
        normalize(v);
    }

	//ogni vertice originale diventa una faccia nel duale
	for (unsigned int i = 0; i < vertices.size(); ++i) {
		vector<int> adjacentFaces;

		// trovo tutte le facce che contengono il vertice i
		for (unsigned int j = 0; j < faces.size(); ++j) {
			if (find(faces[j].vertex_ids.begin(), faces[j].vertex_ids.end(), i) != faces[j].vertex_ids.end()) {
				adjacentFaces.push_back(j);
			}
		}

		if (adjacentFaces.empty()) continue;

		vector<int> orderedFaces;
		set<int> visited;
		visited.insert(adjacentFaces[0]);
		orderedFaces.push_back(adjacentFaces[0]);

		while (true) {
			bool found = false;
			int lastFaceId = orderedFaces.back();
			const face& lastFace = faces[lastFaceId];

			for (int fid : adjacentFaces) {
				if (visited.count(fid)) continue;
				const face& f = faces[fid];

				// controllo se f condivide un lato e il vertice i con lastFace
				int common = 0;
				for (int v1 : f.vertex_ids) {
					for (int v2 : lastFace.vertex_ids) {
						if (v1 == v2) ++common;
					}
				}
				if (common >= 2) { // condividono un lato
					orderedFaces.push_back(fid);
					visited.insert(fid);
					found = true;
					break;
				}
			}

			if (!found) break; // ciclo chiuso o spezzato
		}

		// se il ciclo è incompleto, completo con quelle rimanenti
		for (int fid : adjacentFaces)
			if (!visited.count(fid)) orderedFaces.push_back(fid);

		face df;
		df.id = i;
		for (int fid : orderedFaces) {
			df.vertex_ids.push_back(fid);
		}
		dualFaces.push_back(df);
	}


    // costruzione esplicita degli spigoli nel duale
    map<pair<int, int>, int> edgeMap;
    for (auto& f : dualFaces) {
        int n = f.vertex_ids.size();
        f.edge_ids.clear();
        for (int i = 0; i < n; ++i) {
            int a = f.vertex_ids[i];
            int b = f.vertex_ids[(i + 1) % n];
            auto key = minmax(a, b);
            if (!edgeMap.count(key)) {
                int eid = static_cast<int>(dualEdges.size());
                edgeMap[key] = eid;
                dualEdges.push_back({eid, key.first, key.second});
            }
            f.edge_ids.push_back(edgeMap[key]);
        }
    }

    dualPoly.id = original.id + 1000;
    dualPoly.vertex_ids.clear();
    dualPoly.edge_ids.clear();
    dualPoly.face_ids.clear();

    for (const auto& v : dualVertices)
        dualPoly.vertex_ids.push_back(v.id);
    for (const auto& e : dualEdges)
        dualPoly.edge_ids.push_back(e.id);
    for (const auto& f : dualFaces)
        dualPoly.face_ids.push_back(f.id);

    dualPoly.num_vertices = static_cast<int>(dualVertices.size());
    dualPoly.num_edges = static_cast<int>(dualEdges.size());
    dualPoly.num_faces = static_cast<int>(dualFaces.size());
}

void exportCell0Ds (const vector<vertex>& vertices, const string& filename) { // O(n)
	ofstream file(filename);
	if (!file.is_open()){
		cerr << "Errore nell'apertura del file Cell0Ds" << filename << endl;
		return;
	}
	
	for (const auto &v : vertices){
		file << "ID: " << v.id << "\n";
		file << "Coordinate: " << v.x << " " << v.y << " " << v.z << "\n";
		//file << "ShortPath: " << v.ShortPath << "\n\n";
	}
	file.close();
}

void exportCell1Ds (const vector<edge>& edges, const string& filename)
{
	ofstream file(filename);
	if (!file.is_open()){
		cerr << "Errore nell'apertura del file Cell1Ds" << filename << endl;
		return;
	}
	
	for (const auto &e : edges){
		file << "ID: " << e.id << "\n";
		file << "Origin: " << e.origin << "\n";
		file << "End: " << e.end << "\n";		
	}
	file.close();
}

void exportCell2Ds (const vector<face>& faces, const string& filename)
{
	ofstream file(filename);
	if (!file.is_open()){
		cerr << "Errore nell'apertura del file Cell2Ds" << filename << endl;
		return;
	}
	
	for (const auto &f : faces) {
        file << "ID: " << f.id << "\n";
		file << "Num Vertici: " << f.vertex_ids.size() << "\n";
		file << "Num Lati: " << f.edge_ids.size() << "\n";
		file << "Vertici: ";
		for (auto v : f.vertex_ids) file << v << " ";
		file << "\n";
		file << "Lati: ";
		for (auto e : f.edge_ids) file << e << " ";
		file << "\n";
	}
	file.close();
}

void exportCell3Ds(const vector<polyhedron>& polyhedra, const string& filename)
{
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Errore nell'apertura del file Cell3Ds: " << filename << endl;
        return;
    }

    for (const auto& p : polyhedra) {
        file << "ID: " << p.id << "\n";
        file << "Num Vertici: " << p.vertex_ids.size() << "\n";
        file << "Num Lati: " << p.edge_ids.size() << "\n";
        file << "Num Facce: " << p.face_ids.size() << "\n";

        file << "Vertici: ";
        for (auto v : p.vertex_ids) file << v << " ";
        file << "\n";

        file << "Lati: ";
        for (auto e : p.edge_ids) file << e << " ";
        file << "\n";

        file << "Facce: ";
        for (auto f : p.face_ids) file << f << " ";
        file << "\n\n";
    }

    file.close();
}

void exportToParaview(const vector<vertex>& vertices, const vector<edge>& edges, const string& outputDirectory) { // O(n+m)
    using namespace Gedim;

    UCDUtilities exporter;

    // Vertici
    Eigen::MatrixXd points(3, vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        points(0, i) = vertices[i].x;
        points(1, i) = vertices[i].y;
        points(2, i) = vertices[i].z;
    }

    Eigen::VectorXi pointMaterials = Eigen::VectorXi::Zero(vertices.size());

    UCDProperty<double> vertex_prop;
    vertex_prop.NumComponents = 1;
    vertex_prop.Label = "Visited Nodes";
    vertex_prop.UnitLabel = "";
    vector<double> vdata(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i)
        vdata[i] = vertices[i].ShortPath;
    vertex_prop.Data = vdata.data();
    vector<UCDProperty<double>> pointProps = { vertex_prop };

    string pointFile = outputDirectory + "/vertices.inp";
    exporter.ExportPoints(pointFile, points, pointProps, pointMaterials);
	
	Eigen::MatrixXi segments(2, edges.size());
	for (size_t i = 0; i < edges.size(); ++i) {
		segments(0, i) = edges[i].origin;
		segments(1, i) = edges[i].end;
	}

	Eigen::VectorXi edgeMaterials = Eigen::VectorXi::Zero(edges.size());

	UCDProperty<double> edge_prop;
	edge_prop.NumComponents = 1;
	edge_prop.Label = "Visited Edges";
	edge_prop.UnitLabel = "";
	vector<double> edata(edges.size());
	for (size_t i = 0; i < edges.size(); ++i)
		edata[i] = edges[i].ShortPath;
	edge_prop.Data = edata.data();
	vector<UCDProperty<double>> edgeProps = { edge_prop };

	string edgeFile = outputDirectory + "/edges.inp";
	exporter.ExportSegments(edgeFile, points, segments, pointProps, edgeProps, edgeMaterials);
}

void buildDualFromBaseThenGeodesic(int p, int q, int b, int c, vector<vertex>& finalVertices, vector<edge>& finalEdges, vector<face>& finalFaces, polyhedron& finalPoly) {
	if ((b == 0 || c == 0) && (b + c) > 0 && q != 3) {
		cerr << "Errore.\n";
		return;
	}

	if (p == 3 && q == 3 && b == 0 && c == 0) {
        buildPolyhedron(p, q, b, c, finalVertices, finalEdges, finalFaces, finalPoly);
        return;
    }
    vector<vertex> baseVertices;
    vector<edge> baseEdges;
    vector<face> baseFaces;
    polyhedron basePoly;
    buildPolyhedron(p, q, 0, 0, baseVertices, baseEdges, baseFaces, basePoly);

    vector<vertex> dualVertices;
    vector<edge> dualEdges;
    vector<face> dualFaces;
    polyhedron dualPoly;
    buildDualPolyhedron(baseVertices, baseEdges, baseFaces, basePoly, dualVertices, dualEdges, dualFaces, dualPoly);

	bool classII = (b == c && b > 0);
	bool classI  = !classII && ((b == 0 && c > 0) || (b > 0 && c == 0));


	if (classII) {
		buildClassIIGeodesic(3, p, b, c, dualVertices, dualEdges, dualFaces, dualPoly);
		finalVertices = dualVertices;
		finalEdges    = dualEdges;
		finalFaces    = dualFaces;
		finalPoly     = dualPoly;
	} else if (classI) {
		buildClassIGeodesic(3, p, b + c, dualVertices, dualEdges, dualFaces, dualPoly);
		finalVertices = dualVertices;
		finalEdges    = dualEdges;
		finalFaces    = dualFaces;
		finalPoly     = dualPoly;
    } else if (b == 0 && c == 0) {
        finalVertices = dualVertices;
        finalEdges = dualEdges;
        finalFaces = dualFaces;
        finalPoly = dualPoly;
    } else {
        cerr << "Parametri (b=" << b << ", c=" << c << ") non validi per una geodesica.\n";
        finalVertices = dualVertices;
        finalEdges = dualEdges;
        finalFaces = dualFaces;
        finalPoly = dualPoly;
    }
}