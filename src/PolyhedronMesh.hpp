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
    vertex(): id(-1),x(0.0),y(0.0),z(0.0){}
    vertex(int _id, double _x, double _y, double _z) : id(_id), x(_x), y(_y), z(_z){}
    //operatore somma vettoriale
    vertex operator+(const vertex& other) const {
        return {-1, x + other.x, y + other.y, z + other.z};
    }
    //moltiplicazione per scalare
    vertex operator*(double scalar) const{
        return{-1, x * scalar, y * scalar, z * scalar };
    }
    // sottrazione vettoriale
    vertex operator-(const vertex& other)const {
        return {-1, x - other.x, y - other.y, z - other.z};
    }
    //prodotto scalare
    double dot(const vertex& other) const {
        return  x * other.x + y * other.y + z * other.z;
    }
    // prodotto vettoriale
    vertex cross(const vertex& other) const {
        return {-1, y * other.z - z * other.y, z * other.x - x*other.z, x *other.y - y*other.x};
    }
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
    edge() : id(-1), origin(-1), end(-1) {}
    edge(int _id, int _origin, int _end) : id(_id), origin(_origin), end(_end) {}
};

struct face {
    int id;
    std::vector<int> vertex_ids;
    std::vector<int> edge_ids;
    face() : id(-1) {}
    face(int _id) : id(_id) {}
	
	face(int _id, const std::vector<int>& verts, const std::vector<int>& edges)
        : id(_id), vertex_ids(verts), edge_ids(edges) {}
};

struct polyhedron {
    int id;
    int num_vertices;
    int num_edges;
    int num_faces;
    std::vector<int> vertex_ids;
    std::vector<int> edge_ids;
    std::vector<int> face_ids;
    polyhedron() : id(-1), num_vertices(0), num_edges(0), num_faces(0) {}
    polyhedron(int _id) : id(_id), num_vertices(0), num_edges(0), num_faces(0) {}
};
}
