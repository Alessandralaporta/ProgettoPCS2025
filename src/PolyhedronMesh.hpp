#pragma once
#include <vector>
#include <iostream>
#include <map>
#include <cmath>
#include <algorithm>

namespace PolyhedronMesh {

struct vertex {
    int id;
    double x, y, z;
	int ShortPath = 0;
    //lunghezza del vettore
    double length() const {
        return sqrt(x*x + y*y + z*z);
    }
    // normalizza il vettore e lo trasforma in un versore
    vertex normalized() const {
        double len = length();
        if (len > 1e-12) {
            return {-1, x / len, y / len, z / len};
        }
        return *this;
    }
};

struct edge {
    int id;
    int origin;
    int end;
    
    double length = 0.0;
    int ShortPath = 0;
};

struct face {
    int id;
    std::vector<int> vertex_ids;
    std::vector<int> edge_ids;
};

struct polyhedron {
    int id;
    int num_vertices;
    int num_edges;
    int num_faces;
    std::vector<int> vertex_ids;
    std::vector<int> edge_ids;
    std::vector<int> face_ids;
};
}
