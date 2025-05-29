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


using namespace std;
using namespace PolyhedronMesh;

//ordinare in cell2Ds
//tipo di polyhedron

void normalize(vertex& v) {
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len == 0)  {
		return;
	}	
    v.x /= len;
    v.y /= len;
    v.z /= len;
}

void buildTetrahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
	
	vertices.clear(); edges.clear(); faces.clear();
	
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
		{0, {0, 1, 2}, {0, 1, 3}},
		{1, {0, 1, 3}, {0, 2, 4}},
		{2, {0, 2, 3}, {1, 2, 5}},
		{3, {1, 2, 3}, {3, 4, 5}}
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

void buildEsahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
	
	vertices.clear(); edges.clear(); faces.clear();
	
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
	
	vertices.clear(); edges.clear(); faces.clear();
	
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
        {4, {2, 0, 5}, {0, 8, 9}},  
        {5, {1, 2, 5}, {1, 9,10}},
        {6, {3, 1, 5}, {2,10,11}},
        {7, {0, 3, 5}, {3,11, 8}}
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

const double phi = (1.0 + std::sqrt(5.0)) / 2.0;

void buildDodecahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
	
	vertices.clear(); edges.clear(); faces.clear();
	
    const float a = 1.0 / phi;
    const float b = 1.0;

    vertices = {
        { 0,  b,  b,  b}, { 1,  b,  b, -b}, { 2,  b, -b,  b}, { 3,  b, -b, -b},
        { 4, -b,  b,  b}, { 5, -b,  b, -b}, { 6, -b, -b,  b}, { 7, -b, -b, -b},
        { 8,  0,  a,  phi}, { 9,  0,  a, -phi}, {10,  0, -a,  phi}, {11,  0, -a, -phi},
        {12,  a,  phi, 0}, {13,  a, -phi, 0}, {14, -a,  phi, 0}, {15, -a, -phi, 0},
        {16,  phi, 0,  a}, {17,  phi, 0, -a}, {18, -phi, 0,  a}, {19, -phi, 0, -a}
    };

    edges = {
        { 0, 0, 8}, { 1, 0,12}, { 2, 0, 4}, { 3, 0,16}, { 4, 0, 2},
        { 5, 1, 3}, { 6, 1,12}, { 7, 1, 9}, { 8, 1,17}, { 9, 1,16},
        {10, 2, 3}, {11, 2,13}, {12, 2,10}, {13, 2,16},
        {14, 3,11}, {15, 3,17},
        {16, 4, 8}, {17, 4,14}, {18, 4,18}, {19, 5, 9}, {20, 5,14},
        {21, 5,19}, {22, 6,10}, {23, 6,15}, {24, 6,18}, {25, 7,11},
        {26, 7,15}, {27, 7,19}, {28, 8,10}, {29, 9,11}
    };

    faces = {
        {0,  {0, 4,13,12,28},  {0, 4,13,12,28}}, 
        {1,  {1, 6,7,8,9},     {1, 6,7,8,9}},     
        {2,  {2,16,17,18,1},   {2,16,17,18,1}},
        {3,  {3,9,13,11,10},   {3,9,13,11,10}},
        {4,  {5,10,15,14,6},   {5,10,15,14,6}},
        {5,  {7,19,20,17,16},  {7,19,20,17,16}},
        {6,  {8,19,29,25,14},  {8,19,29,25,14}},
        {7,  {11,12,22,23,10}, {11,12,22,23,10}},
        {8,  {15,26,27,21,20}, {15,26,27,21,20}},
        {9,  {18,24,23,22,16}, {18,24,23,22,16}},
        {10, {5,3,13,11,14},   {5,3,13,11,14}}, 
        {11, {27,26,25,24,28}, {27,26,25,24,28}} 
    };
	
	polyhedron.id = 3;
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

void buildIcosahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
	
	vertices.clear(); edges.clear(); faces.clear();
	
	 vertices = {
        {0, -1,  phi, 0}, {1, 1,  phi, 0}, {2, -1, -phi, 0}, {3, 1, -phi, 0},
        {4, 0, -1,  phi}, {5, 0, 1,  phi}, {6, 0, -1, -phi}, {7, 0, 1, -phi},
        {8,  phi, 0, -1}, {9,  phi, 0, 1}, {10, -phi, 0, -1}, {11, -phi, 0, 1}
    };

    edges = {
        {0, 0, 1}, {1, 0, 5}, {2, 0, 11}, {3, 0, 4}, {4, 0, 10},
        {5, 1, 5}, {6, 1, 9}, {7, 1, 8}, {8, 1, 4}, {9, 2, 3},
        {10, 2, 7}, {11, 2, 6}, {12, 2, 4}, {13, 2, 10}, {14, 3, 6},
        {15, 3, 8}, {16, 3, 7}, {17, 3, 9}, {18, 4, 11}, {19, 5, 11},
        {20, 6, 10}, {21, 6, 8}, {22, 7, 6}, {23, 7, 10}, {24, 8, 9},
        {25, 9, 5}, {26, 11, 10}, {27, 7, 5}, {28, 8, 6}, {29, 9, 3}
    };

    faces = {
        {0, {0, 1, 5}, {0, 1, 5}}, {1, {0, 4, 3}, {0, 4, 3}}, {2, {0, 3, 1}, {0, 3, 1}},
        {3, {1, 3, 5}, {1, 3, 5}}, {4, {2, 4, 0}, {2, 4, 0}}, {5, {2, 0, 1}, {2, 0, 1}},
        {6, {3, 0, 4}, {3, 0, 4}}, {7, {3, 1, 5}, {3, 1, 5}}, {8, {4, 0, 1}, {4, 0, 1}},
        {9, {5, 1, 3}, {5, 1, 3}}, {10, {6, 2, 10}, {11, 13, 20}}, {11, {6, 10, 3}, {20, 14, 15}},
        {12, {6, 3, 8}, {14, 15, 21}}, {13, {7, 2, 6}, {10, 11, 22}}, {14, {7, 6, 8}, {22, 21, 28}},
        {15, {8, 3, 9}, {15, 17, 24}}, {16, {8, 9, 5}, {24, 25, 5}}, {17, {9, 3, 1}, {17, 7, 6}},
        {18, {9, 1, 5}, {6, 5, 25}}, {19, {10, 2, 4}, {13, 12, 18}}
    };

    polyhedron.id = 4;
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

void buildPolyhedron(int p, int q, int b, int c, std::vector<vertex> &vertices, std::vector<edge> &edges, std::vector<face> &faces, polyhedron &polyhedron) {
	if (p < 3 || q < 3) {
		cout << "Valori non validi" << endl;
	}
	
	bool ClassI = (b == 0 && c > 0) || (b > 0 && c == 0);
	bool ClassII = (b == c && b > 0);
	
	if (!ClassI && !ClassII && (b != 0 || c != 0)) {
		cout << "Tipo di classe non supportato" << endl;
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
}

bool sameVertex(const vertex& a, const vertex& b, double tolerance) {
    return fabs(a.x - b.x) < tolerance && fabs(a.y - b.y) < tolerance && fabs(a.z - b.z) < tolerance;
}


int getOrAddVertex(double x, double y, double z, vector<vertex>& vertices) {
    for (auto& existing : vertices) {
        if (sameVertex(existing, vertex{-1, x, y, z})) {
            // Se il vertice è già presente MA non ha ancora un id valido, lo assegna ora
            if (existing.id == -1)
                existing.id = &existing - &vertices[0]; // o usare un contatore globale
            return existing.id;
        }
    }
    int newId = vertices.size();
    vertex v{newId, x, y, z};
    vertices.push_back(v);
    return newId;
}


void projectVerticesOnUnitSphere(vector<vertex>&vertices){
    for(auto& v : vertices){
        normalize(v);
    }
}

void buildGeodesicPolyhedron(int p, int q, int b, int c, vector<vertex>& vertices, vector<edge>& edges, vector<face>& faces, polyhedron& poly) {
    if ((p < 3 || p > 5) || (q != 3 && q != 4 && q != 5)) {
        cerr << "Tipo non supportato: p = " << p << ", q = " << q << endl;
        return;
    }

    vector<vertex> baseVertices;
    vector<edge> baseEdges;
    vector<face> baseFaces;
    polyhedron basePoly;

    buildPolyhedron(p, q, 0, 0, baseVertices, baseEdges, baseFaces, basePoly);

    // triangola tutte le facce (se necessario)
    vector<face> triangleFaces;
    int faceId = 0;
    for (const auto& f : baseFaces) {
        const auto& v = f.vertex_ids;
        if (v.size() == 3) {
            triangleFaces.push_back({faceId++, v, {}});
        } else if (v.size() > 3) {
            for (size_t i = 1; i < v.size() - 1; ++i) {
                triangleFaces.push_back({faceId++, {v[0], v[i], v[i + 1]}, {}});
            }
        } else {
           cerr << "Faccia con meno di 3 vertici, ignorata." << endl;
        }
    }

    // suddivisione geodetica
    int N = b + c;
    for (const auto& f : triangleFaces) {
        vertex A = baseVertices[f.vertex_ids[0]];
        vertex B = baseVertices[f.vertex_ids[1]];
        vertex C = baseVertices[f.vertex_ids[2]];

        vector<vector<int>> grid(N + 1);

        for (int i = 0; i <= N; ++i) {
            grid[i].resize(i + 1);
            for (int j = 0; j <= i; ++j) {
                double alpha = (double)(b - j) / N;
                double beta = (double)(c - (i - j)) / N;
                double gamma = 1.0 - alpha - beta;

                double x = alpha * A.x + beta * B.x + gamma * C.x;
                double y = alpha * A.y + beta * B.y + gamma * C.y;
                double z = alpha * A.z + beta * B.z + gamma * C.z;

                int id = getOrAddVertex(x, y, z, vertices);
                grid[i][j] = id;
            }
        }

        // crea triangoli dalla griglia
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < i + 1; ++j) {
                int v0 = grid[i][j];
                int v1 = grid[i + 1][j];
                int v2 = grid[i + 1][j + 1];
                faces.push_back({(int)faces.size(), {v0, v1, v2}, {}});

                if (j < i) {
                    int v3 = grid[i][j + 1];
                    faces.push_back({(int)faces.size(), {v0, v2, v3}, {}});
                }
            }
        }
    }

    // crea gli spigoli
    edges.clear();
    map<int, int> edgeMap;
    const int MAX_VERTS = 10000;

    for (auto& f : faces) {
        for (int i = 0; i < 3; ++i) {
            int a = f.vertex_ids[i];
            int b = f.vertex_ids[(i + 1) % 3];
            int key = min(a, b) * MAX_VERTS + std::max(a, b);

            if (edgeMap.count(key) == 0) {
                int id = edges.size();
                edges.push_back({id, a, b});
                edgeMap[key] = id;
            }
            f.edge_ids.push_back(edgeMap[key]);
        }
    }
	
    poly.vertex_ids.clear();
    for (const auto& v : vertices) poly.vertex_ids.push_back(v.id);
    poly.edge_ids.clear();
    for (const auto& e : edges) poly.edge_ids.push_back(e.id);
    poly.face_ids.clear();
    for (const auto& f : faces) poly.face_ids.push_back(f.id);
	
	poly.id = basePoly.id;
}
	
bool isFaceConsistent(const face& face, const std::vector<edge>& edges) {
    const auto& ids = face.edge_ids;
    if (ids.size() < 3) return false;

    for (size_t i = 0; i < ids.size(); ++i) {
        size_t current = ids[i];
        size_t next = ids[(i + 1) % ids.size()];
        
        if (current >= edges.size() || next >= edges.size())
            return false;

        if (edges[current].end != edges[next].origin)
            return false; 
    }
    return true;
}

double distance(const vertex& a, const vertex &b) {
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z + b.z));
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

    // Costruzione grafo come lista adiacenza
    vector<vector<pair<int, double>>> adj(N);
    for (const auto& e : edges) {
        double len = e.length > 0 ? e.length : distance(vertices[e.origin], vertices[e.end]);
        adj[e.origin].emplace_back(e.end, len);
        adj[e.end].emplace_back(e.origin, len); // non orientato
    }

    // Algoritmo di Dijkstra
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

    // Ricostruzione cammino
    vector<int> path;
    for (int at = endId; at != -1; at = prev[at])
        path.push_back(at);
    reverse(path.begin(), path.end());

    // Resetta ShortPath
    for (auto& v : vertices) v.ShortPath = 0;
    for (auto& e : edges) e.ShortPath = 0;

    // Marca i vertici
    for (int v : path)
        vertices[v].ShortPath = 1;

    // Marca i lati e somma lunghezza
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
}

void exportCell0Ds (const vector<vertex>& vertices, const string& filename) {
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

void buildDualPolyhedron(const vector<vertex>& vertices, const vector<face>& faces, const polyhedron& original, vector<vertex>& dualVertices, vector<face>& dualFaces, polyhedron& dualPoly) {

    dualVertices = calculateCentroids(vertices, faces);
    map<int, vector<int>> vertexToFaceIds;

    // Mappa ogni vertice alle facce che lo contengono
    for (const auto& f : faces) {
        for (int vid : f.vertex_ids) {
            vertexToFaceIds[vid].push_back(f.id);
        }
    }

    int newFaceId = 0;
    dualFaces.clear();
    dualFaces.reserve(vertices.size());

    for (const auto& [vId, adjacentFaceIds] : vertexToFaceIds) {
        const vertex& center = vertices[vId];
        vector<pair<double, int>> ordered;

        for (int fid : adjacentFaceIds) {
            const vertex& fcenter = dualVertices[fid];
            double dx = fcenter.x - center.x;
            double dy = fcenter.y - center.y;
            double angle = atan2(dy, dx);
            ordered.emplace_back(angle, fid);
        }

        sort(ordered.begin(), ordered.end());

        face newFace(newFaceId++);
        for (const auto& [angle, fid] : ordered)
            newFace.vertex_ids.push_back(fid);

        dualFaces.push_back(newFace);
    }

    dualPoly.id = original.id + 1000;
    dualPoly.vertex_ids.clear();
    dualPoly.edge_ids.clear();
    dualPoly.face_ids.clear();

    for (const auto& v : dualVertices)
        dualPoly.vertex_ids.push_back(v.id);
    for (const auto& f : dualFaces)
        dualPoly.face_ids.push_back(f.id);

    // Costruzione degli spigoli unici del duale
    set<pair<int, int>> uniqueEdges;
    int edgeId = 0;

    for (const auto& f : dualFaces) {
        int n = f.vertex_ids.size();
        for (int i = 0; i < n; ++i) {
            int a = f.vertex_ids[i];
            int b = f.vertex_ids[(i + 1) % n];

            if (a > b) swap(a, b);
            if (uniqueEdges.insert({a, b}).second) {
                dualPoly.edge_ids.push_back(edgeId++);
            }
        }
    }

    dualPoly.num_vertices = dualVertices.size();
    dualPoly.num_faces = dualFaces.size();
    dualPoly.num_edges = dualPoly.edge_ids.size();
}

