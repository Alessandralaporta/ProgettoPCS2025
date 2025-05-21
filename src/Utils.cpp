#include "Utils.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "PolyhedronMesh.hpp"

using namespace std;
using namespace PolyhedronMesh;

//ordinare in cell2Ds
//tipo di polyhedron



void normalize(vertex &v) {
    double norm = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	if (norm > 1e-12) {  
        v.x /= norm;
        v.y /= norm;
        v.z /= norm;
    }
}

void buildTetrahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
	
	//vertices.clear(); edges.clear(); faces.clear();
	
	vertices = {
		{0, 1, 1, 1},
		{1, -1, -1, 1},
		{2, -1, 1, -1},
		{3, 1, -1, -1}
	};
	
	for (size_t i = 0; i < vertices.size(); ++i) {
		normalize(vertices[i]);
	}
	
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
	//polyhedron.vertex_ids.clear();
	for (size_t i = 0; i < vertices.size(); ++i) {
		polyhedron.vertex_ids.push_back(i);
	}

	//polyhedron.edge_ids.clear();
	for (size_t i = 0; i < edges.size(); ++i) {
		polyhedron.edge_ids.push_back(i);
	}

	//polyhedron.face_ids.clear();
	for (size_t i = 0; i < faces.size(); ++i) {
		polyhedron.face_ids.push_back(i);
	}
}

void buildEsahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
	
	//vertices.clear(); edges.clear(); faces.clear();
	
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
	
	for (size_t i = 0; i < vertices.size(); ++i) {
		normalize(vertices[i]);
	}
	
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
    //polyhedron.vertex_ids.clear();
	for (size_t i = 0; i < vertices.size(); ++i) {
		polyhedron.vertex_ids.push_back(i);
	}

	//polyhedron.edge_ids.clear();
	for (size_t i = 0; i < edges.size(); ++i) {
		polyhedron.edge_ids.push_back(i);
	}

	//polyhedron.face_ids.clear();
	for (size_t i = 0; i < faces.size(); ++i) {
		polyhedron.face_ids.push_back(i);
	}
}

void buildOctahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
	
	//vertices.clear(); edges.clear(); faces.clear();
	
	 vertices = {
        {0,  1,  0,  0},
        {1, -1,  0,  0},  
        {2,  0,  1,  0},  
        {3,  0, -1,  0},  
        {4,  0,  0,  1}, 
        {5,  0,  0, -1}   
    };
	
	for (size_t i = 0; i < vertices.size(); ++i) {
		normalize(vertices[i]);
	}
	
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
    //polyhedron.vertex_ids.clear();
	for (size_t i = 0; i < vertices.size(); ++i) {
		polyhedron.vertex_ids.push_back(i);
	}

	//polyhedron.edge_ids.clear();
	for (size_t i = 0; i < edges.size(); ++i) {
		polyhedron.edge_ids.push_back(i);
	}

	//polyhedron.face_ids.clear();
	for (size_t i = 0; i < faces.size(); ++i) {
		polyhedron.face_ids.push_back(i);
	}
	
}

const double phi = (1.0 + std::sqrt(5.0)) / 2.0;

void buildDodecahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
	
	//vertices.clear(); edges.clear(); faces.clear();
	
    const float a = 1.0 / phi;
    const float b = 1.0;

    vertices = {
        { 0,  b,  b,  b}, { 1,  b,  b, -b}, { 2,  b, -b,  b}, { 3,  b, -b, -b},
        { 4, -b,  b,  b}, { 5, -b,  b, -b}, { 6, -b, -b,  b}, { 7, -b, -b, -b},
        { 8,  0,  a,  phi}, { 9,  0,  a, -phi}, {10,  0, -a,  phi}, {11,  0, -a, -phi},
        {12,  a,  phi, 0}, {13,  a, -phi, 0}, {14, -a,  phi, 0}, {15, -a, -phi, 0},
        {16,  phi, 0,  a}, {17,  phi, 0, -a}, {18, -phi, 0,  a}, {19, -phi, 0, -a}
    };

    for (size_t i = 0; i < vertices.size(); ++i) {
        normalize(vertices[i]);
    }

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
	//polyhedron.vertex_ids.clear();
	for (size_t i = 0; i < vertices.size(); ++i) {
		polyhedron.vertex_ids.push_back(i);
	}

	//polyhedron.edge_ids.clear();
	for (size_t i = 0; i < edges.size(); ++i) {
		polyhedron.edge_ids.push_back(i);
	}

	//polyhedron.face_ids.clear();
	for (size_t i = 0; i < faces.size(); ++i) {
		polyhedron.face_ids.push_back(i);
	}
}

void buildIcosahedron(vector<vertex> &vertices, vector<edge> &edges, vector<face> &faces, polyhedron &polyhedron) {
	
	//vertices.clear(); edges.clear(); faces.clear();
	
	 vertices = {
        {0, -1,  phi, 0}, {1, 1,  phi, 0}, {2, -1, -phi, 0}, {3, 1, -phi, 0},
        {4, 0, -1,  phi}, {5, 0, 1,  phi}, {6, 0, -1, -phi}, {7, 0, 1, -phi},
        {8,  phi, 0, -1}, {9,  phi, 0, 1}, {10, -phi, 0, -1}, {11, -phi, 0, 1}
    };

    for (size_t i = 0; i < vertices.size(); ++i) {
        normalize(vertices[i]);
    }

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
    //polyhedron.vertex_ids.clear();
	for (size_t i = 0; i < vertices.size(); ++i) {
		polyhedron.vertex_ids.push_back(i);
	}

	//polyhedron.edge_ids.clear();
	for (size_t i = 0; i < edges.size(); ++i) {
		polyhedron.edge_ids.push_back(i);
	}

	//polyhedron.face_ids.clear();
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
	
bool isFaceConsistent(const face& face, const std::vector<edge>& edges) {
    for (int i = 0; i < face.edge_ids.size(); ++i) {
        int current = face.edge_ids[i];
        int next = face.edge_ids[(i + 1) % face.edge_ids.size()];
        if (edges[current].end != edges[next].origin)
            return false; 
    }
    return true;

}

void exportCell0Ds (const vector<vertex>& vertices, const string& filename) {
	ofstream file(filename);
	if (!file.is_open()){
		cerr << "Errore nell'apertura del file Cell0Ds" << filename << endl;
		return;
	}
	
	for (const auto &v : vertices){
		file << v.id << " " << v.x << " " << v.y << " " << v.z << " " << v.ShortPath << endl;		
	}
	file.close();
}

void exportCell1Ds (const vector<edge>& edges, const string& filename = "cell1Ds.txt")
{
	ofstream file(filename);
	if (!file.is_open()){
		cerr << "Errore nell'apertura del file Cell1Ds" << filename << endl;
		return;
	}
	
	for (const auto &e : edges){
		file << e.id << " " << e.origin << " " << e.end << endl;		
	}
	file.close();
}

void exportCell2Ds (const vector<face>& faces, const string& filename = "cell2Ds.txt")
{
	ofstream file(filename);
	if (!file.is_open()){
		cerr << "Errore nell'apertura del file Cell2Ds" << filename << endl;
		return;
	}
	
	for (const auto &f : faces) {
        file << f.id << " " << f.vertex_ids.size() << " " << f.edge_ids.size() << " ";
        for (auto v : f.vertex_ids) file << v << " ";
        for (auto e : f.edge_ids) file << e << " ";
        file << endl;
	}
	file.close();
}

void exportCell3Ds (const vector<polyhedron>& polyhedra, const string& filename = "cell3Ds.txt")
{
	ofstream file(filename);
	if (!file.is_open()){
		cerr << "Errore nell'apertura del file Cell3Ds" << filename << endl;
		return;
	}
	
	for (const auto &p : polyhedra) {
        file << p.id << " " << p.vertex_ids.size() << " " << p.edge_ids.size() << " " << p.face_ids.size() << " ";
        for (auto v : p.vertex_ids) file << v << " ";
        for (auto e : p.edge_ids) file << e << " ";
        for (auto f : p.face_ids) file << f << " ";
        file << endl;
    }
	file.close();
}
