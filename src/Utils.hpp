#pragma once
#include <string>
#include "PolyhedronMesh.hpp"

void exportCell0Ds (const vector<vertex>& vertices, const string& filename = "cell0Ds.txt");

void exportCell1Ds (const vector<edge>& edges, const string& filename = "cell1Ds.txt");

void exportCell2Ds (const vector<face>& faces, const string& filename = "cell2Ds.txt");

void exportCell3Ds (const vector<polyhedron>& polyhedra, const string& filename = "cell3Ds.txt");


