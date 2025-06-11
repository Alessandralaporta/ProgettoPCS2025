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

void normalize(vertex& v) {
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len < 1e-8) return; //tolleranza
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
		{0, {0, 1, 2}, {0, 3, 1}},   // edge 0: 0→1, edge 3: 1→2, edge 1: 2→0
		{1, {0, 2, 3}, {1, 5, 2}},   // edge 1: 0→2, edge 5: 2→3, edge 2: 3→0
		{2, {0, 3, 1}, {2, 4, 0}},   // edge 2: 0→3, edge 4: 3→1, edge 0: 1→0
		{3, {1, 2, 3}, {3, 5, 4}}    // edge 3: 1→2, edge 5: 2→3, edge 4: 3→1
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
		{0, {0, 1, 2, 3},  {0, 1, 2, 3}},     // base inferiore
		{1, {4, 5, 6, 7},  {4, 5, 6, 7}},     // base superiore
		{2, {0, 1, 5, 4},  {0, 9, 4, 8}},     // lato fronte
		{3, {1, 2, 6, 5},  {1,10, 5, 9}},     // lato destra
		{4, {2, 3, 7, 6},  {2,11, 6,10}},     // lato retro
		{5, {3, 0, 4, 7},  {3, 8, 7,11}}      // lato sinistra
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
		{0, {0, 2, 4}, {0, 5, 4}},     // 0→2, 2→4, 4→0
		{1, {2, 1, 4}, {1, 6, 5}},     // 2→1, 1→4, 4→2
		{2, {1, 3, 4}, {2, 7, 6}},     // 1→3, 3→4, 4→1
		{3, {3, 0, 4}, {3, 4, 7}},     // 3→0, 0→4, 4→3
		{4, {0, 5, 2}, {8, 9, 0}},     // 0→5, 5→2, 2→0
		{5, {2, 5, 1}, {9,10, 1}},     // 2→5, 5→1, 1→2
		{6, {1, 5, 3}, {10,11, 2}},    // 1→5, 5→3, 3→1
		{7, {3, 5, 0}, {11, 8, 3}}     // 3→5, 5→0, 0→3
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

const double phi = (1.0 + sqrt(5.0)) / 2.0;

void buildDodecahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
    vertices.clear();
    edges.clear();
    faces.clear();

    const double phi = (1.0 + sqrt(5.0)) / 2.0;
    const double a = 1.0 / phi;
    const double b = 1.0;

    // Vertici numerati da 0 a 19
    vertices = {
        { 0,  b,  b,  b}, { 1,  b,  b, -b}, { 2,  b, -b,  b}, { 3,  b, -b, -b},
        { 4, -b,  b,  b}, { 5, -b,  b, -b}, { 6, -b, -b,  b}, { 7, -b, -b, -b},
        { 8,  0,  a,  phi}, { 9,  0,  a, -phi}, {10,  0, -a,  phi}, {11,  0, -a, -phi},
        {12,  a,  phi, 0}, {13,  a, -phi, 0}, {14, -a,  phi, 0}, {15, -a, -phi, 0},
        {16,  phi, 0,  a}, {17,  phi, 0, -a}, {18, -phi, 0,  a}, {19, -phi, 0, -a}
    };

    // Facce: ogni faccia ha 5 vertici
	vector<vector<int>> vertex_id_faces = {
		{0, 1, 2, 3, 4},
		{0, 5, 10, 6, 1},
		{1, 6, 11, 7, 2},
		{2, 7, 12, 8, 3},
		{3, 8, 13, 9, 4},
		{4, 9, 14, 5, 0},
		{15, 10, 5, 14, 19},
		{16, 11, 6, 10, 15},
		{17, 12, 7, 11, 16},
		{18, 13, 8, 12, 17},
		{19, 14, 9, 13, 18},
		{19, 18, 17, 16, 15}
	};

    // Costruzione delle facce
    for (int i = 0; i < vertex_id_faces.size(); ++i) {
        face f;
        f.id = i;
        f.vertex_ids = vertex_id_faces[i];
        faces.push_back(f);
    }

    // Set per evitare spigoli duplicati
    set<pair<int, int>> edgeSet;

    // Costruzione degli spigoli unici
    for (auto& f : faces) {
        f.edge_ids.clear();
        int n = f.vertex_ids.size();
        for (int i = 0; i < n; ++i) {
            int a = f.vertex_ids[i];
            int b = f.vertex_ids[(i + 1) % n];
            auto key = minmax(a, b);  // Ordina i vertici per evitare duplicati

            // Se l'edge non è stato visto, aggiungilo al set
            if (edgeSet.count(key) == 0) {
                edgeSet.insert(key);  // Aggiungi lo spigolo al set
                edges.push_back({int(edges.size()), key.first, key.second});  // Crea un edge e aggiungilo a 'edges'
            }

            // Trova l'ID dello spigolo corrispondente
            for (size_t j = 0; j < edges.size(); ++j) {
                if (minmax(edges[j].origin, edges[j].end) == key) {
                    f.edge_ids.push_back(j);
                    break;
                }
            }
        }
    }

    // Popolamento finale del polyhedron
    polyhedron.id = 3;
    polyhedron.vertex_ids.clear();
    for (size_t i = 0; i < vertices.size(); ++i) polyhedron.vertex_ids.push_back(i);
    polyhedron.edge_ids.clear();
    for (size_t i = 0; i < edges.size(); ++i) polyhedron.edge_ids.push_back(i);
    polyhedron.face_ids.clear();
    for (size_t i = 0; i < faces.size(); ++i) polyhedron.face_ids.push_back(i);
}

void buildIcosahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
    vertices.clear(); edges.clear(); faces.clear();

    const double phi = (1.0 + sqrt(5.0)) / 2.0;

    // Vertici numerati da 0 a 11
    vertices = {
        { 0, -1,  phi, 0}, { 1,  1,  phi, 0}, { 2, -1, -phi, 0}, { 3,  1, -phi, 0},
        { 4,  0, -1,  phi}, { 5,  0,  1,  phi}, { 6,  0, -1, -phi}, { 7,  0,  1, -phi},
        { 8,  phi, 0, -1}, { 9,  phi, 0,  1}, {10, -phi, 0, -1}, {11, -phi, 0,  1}
    };

    // Definizione delle 20 facce triangolari (solo vertex_ids)
    vector<vector<int>> vertex_id_faces = {
        {0, 1, 5}, {0, 5,11}, {0,11,4}, {0,4,10}, {0,10,1},
        {1,9,5},  {5,9,7},  {5,7,11}, {11,7,3}, {11,3,4},
        {4,3,6},  {4,6,10}, {10,6,2}, {10,2,1}, {1,2,9},
        {9,2,8},  {9,8,7},  {7,8,3},  {3,8,6},  {6,8,2}
    };

    // Costruzione delle facce
    for (int i = 0; i < vertex_id_faces.size(); ++i) {
        face f;
        f.id = i;
        f.vertex_ids = vertex_id_faces[i];
        faces.push_back(f);
    }

    // Set per evitare spigoli duplicati
    set<pair<int, int>> edgeSet;

    // Costruzione degli spigoli unici e assegnazione a ogni faccia
    for (auto& f : faces) {
        f.edge_ids.clear();
        int n = f.vertex_ids.size();
        for (int i = 0; i < n; ++i) {
            int a = f.vertex_ids[i];
            int b = f.vertex_ids[(i + 1) % n];
            auto key = minmax(a, b);  // Ordina per evitare duplicati

            // Aggiungi solo se lo spigolo non esiste ancora
            if (edgeSet.count(key) == 0) {
                int eid = edges.size();
                edges.push_back({eid, key.first, key.second});
                edgeSet.insert(key);
            }

            // Trova e assegna l'ID dello spigolo a questa faccia
            for (const auto& e : edges) {
                if (minmax(e.origin, e.end) == key) {
                    f.edge_ids.push_back(e.id);
                    break;
                }
            }
        }
    }

    // Popolamento finale del polyhedron
    polyhedron.id = 4;
    polyhedron.vertex_ids.clear();
    for (size_t i = 0; i < vertices.size(); ++i) polyhedron.vertex_ids.push_back(i);
    polyhedron.edge_ids.clear();
    for (size_t i = 0; i < edges.size(); ++i) polyhedron.edge_ids.push_back(i);
    polyhedron.face_ids.clear();
    for (size_t i = 0; i < faces.size(); ++i) polyhedron.face_ids.push_back(i);
}

void buildPolyhedron(int p, int q, int b, int c, vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
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
    return distance(a, b) <= tolerance;
}

int getOrAddVertex(double x, double y, double z, std::vector<vertex>& vertices, double tolerance) {
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

void projectVerticesOnUnitSphere(vector<vertex>&vertices){
    for(auto& v : vertices){
        normalize(v);
    }
}

void buildClassIGeodesic(int p, int q, int b, vector<vertex>& vertices, vector<edge>& edges, vector<face>& faces, polyhedron& poly) {
    if ((p < 3 || p > 5) || (q != 3 && q != 4 && q != 5)) {
        cerr << "Tipo non supportato: p = " << p << ", q = " << q << endl;
        return;
    }

    vertices.clear();
    edges.clear();
    faces.clear();

    vector<vertex> baseVertices;
    vector<edge> baseEdges;
    vector<face> baseFaces;
    polyhedron basePoly;

    buildPolyhedron(p, q, 0, 0, baseVertices, baseEdges, baseFaces, basePoly);

    vector<face> triangleFaces;
    int faceId = 0;
    for (const auto& f : baseFaces) {
        if (f.vertex_ids.size() == 3) {
            triangleFaces.push_back({faceId++, f.vertex_ids, {}});
        } else {
            for (size_t i = 1; i < f.vertex_ids.size() - 1; ++i) {
                triangleFaces.push_back({faceId++, {f.vertex_ids[0], f.vertex_ids[i], f.vertex_ids[i + 1]}, {}});
            }
        }
    }

    const int N = b;
    const double tolerance = 1e-5;

    for (const auto& f : triangleFaces) {
        vertex A = baseVertices[f.vertex_ids[0]];
        vertex B = baseVertices[f.vertex_ids[1]];
        vertex C = baseVertices[f.vertex_ids[2]];

        vector<vector<int>> grid(N + 1, vector<int>(N + 1, -1));

        // Crea vertici baricentrici
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
                grid[i][j] = getOrAddVertex(v.x, v.y, v.z, vertices, tolerance);
            }
        }

        // Triangola la griglia baricentrica
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N - i; ++j) {
                int v0 = grid[i][j];
                int v1 = grid[i + 1][j];
                int v2 = grid[i][j + 1];
                if (v0 >= 0 && v1 >= 0 && v2 >= 0)
                    faces.push_back({(int)faces.size(), {v0, v1, v2}, {}});

                if (i + j < N - 1) {
                    int v3 = grid[i + 1][j + 1];
                    if (v1 >= 0 && v2 >= 0 && v3 >= 0)
                        faces.push_back({(int)faces.size(), {v1, v3, v2}, {}});
                }
            }
        }
    }

    // Crea spigoli unici e associa agli ID facce
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

    // Assemble il poly finale
    poly = basePoly;
    poly.vertex_ids.clear();
    for (const auto& v : vertices) poly.vertex_ids.push_back(v.id);
    poly.edge_ids.clear();
    for (const auto& e : edges) poly.edge_ids.push_back(e.id);
    poly.face_ids.clear();
    for (const auto& f : faces) poly.face_ids.push_back(f.id);
    poly.id = basePoly.id;
}


void buildClassIIGeodesic(int p, int q, int b, int c, vector<vertex>& vertices, vector<edge>& edges, vector<face>& faces, polyhedron& poly) {
    if (b <= 0 || c <= 0 || b != c) {
        cerr << "Classe II richiede b = c > 0\n";
        return;
    }

    int N = b + c;
    const double tolerance = 1e-6;
    vector<vertex> baseVertices;
    vector<edge> baseEdges;
    vector<face> baseFaces;
    polyhedron basePoly;
    buildPolyhedron(p, q, 0, 0, baseVertices, baseEdges, baseFaces, basePoly);

    vector<face> triangleFaces;
    int triangleCount = 0;
    for (const auto& f : baseFaces) {
        if (f.vertex_ids.size() < 3) continue;
        for (size_t i = 1; i + 1 < f.vertex_ids.size(); ++i) {
            triangleFaces.push_back({triangleCount++, {f.vertex_ids[0], f.vertex_ids[i], f.vertex_ids[i + 1]}, {}});
        }
    }

    cout << ">> Generati triangoli: " << triangleFaces.size() << endl;

    vertices.clear();
    edges.clear();
    faces.clear();

    for (size_t t = 0; t < triangleFaces.size(); ++t) {
        const auto& f = triangleFaces[t];
        vertex A = baseVertices[f.vertex_ids[0]];
        vertex B = baseVertices[f.vertex_ids[1]];
        vertex C = baseVertices[f.vertex_ids[2]];

        vector<vector<int>> grid(N + 1, vector<int>(N + 1, -1));

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

                int existingId = -1;
                for (size_t vi = 0; vi < vertices.size(); ++vi) {
                    if (sameVertex(vertices[vi], v, tolerance)) {
                        existingId = static_cast<int>(vi);
                        break;
                    }
                }

                if (existingId != -1) {
                    grid[i][j] = existingId;
                } else {
                    int newId = static_cast<int>(vertices.size());
                    v.id = newId;
                    vertices.push_back(v);
                    grid[i][j] = newId;
                }
            }
        }

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N - i; ++j) {
                int v0 = grid[i][j];
                int v1 = grid[i + 1][j];
                int v2 = grid[i][j + 1];
                faces.push_back({static_cast<int>(faces.size()), {v0, v1, v2}, {}});

                if (i + j < N - 1) {
                    int v3 = grid[i + 1][j + 1];
                    faces.push_back({static_cast<int>(faces.size()), {v1, v3, v2}, {}});
                }
            }
        }
    }

    for (auto& f : faces) {
        f.edge_ids.clear();
        for (int i = 0; i < 3; ++i) {
            int a = f.vertex_ids[i];
            int b = f.vertex_ids[(i + 1) % 3];
            bool found = false;
            for (auto& e : edges) {
                if ((e.origin == a && e.end == b) || (e.origin == b && e.end == a)) {
                    f.edge_ids.push_back(e.id);
                    found = true;
                    break;
                }
            }
            if (!found) {
                int eid = static_cast<int>(edges.size());
                edges.push_back({eid, a, b});
                f.edge_ids.push_back(eid);
            }
        }
    }

    poly = basePoly;
    poly.vertex_ids.clear();
    for (const auto& v : vertices) poly.vertex_ids.push_back(v.id);
    poly.edge_ids.clear();
    for (const auto& e : edges) poly.edge_ids.push_back(e.id);
    poly.face_ids.clear();
    for (const auto& f : faces) poly.face_ids.push_back(f.id);
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

double distance(const vertex& a, const vertex& b) {
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
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
	
	for (const auto& v : vertices) {
		if (vertexToFaceIds.find(v.id) == vertexToFaceIds.end()) {
			cerr << "❗Vertice " << v.id << " non è presente in nessuna faccia!\n";
		}
	}
	
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

void exportToParaview(const vector<vertex>& vertices, const vector<edge>& edges, const string& outputDirectory) {
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