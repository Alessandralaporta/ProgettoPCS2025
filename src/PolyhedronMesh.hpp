#pragma once
#include <vector>
#include <Eigen/Eigen>
#include <iostream>

using namespace std;
using namespace Eigen;

namespace PolyhedronMesh{
struct vertex {
	int id;
	double x, y, z;
	int ShortPath =0;				//controllare se int o bool
	void normalize() 
	{
        double length = sqrt(x*x + y*y + z*z);
        x /= length; y /= length; z /= length;
    };
}

struct edge {
	int id;
	int origin;
	int end;
	double length = 0.0;
	int ShortPath = 0;				//controllare se int o bool
}

struct face {
	int id;
	vector<int> vertex_ids;
	vector<int> edge_ids;
}

struct polyhedron {
	int id;
	vector<int> vertex_ids;
	vector<int> edge_ids;
	vector>int> face_ids;
}

}