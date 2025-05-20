#pragma once
#include <string>
#include <vector>
#include "PolyhedronMesh.hpp"

void exportCell0Ds(const std::vector<vertex>& vertices, const std::string& filename = "cell0Ds.txt");

void exportCell1Ds(const std::vector<edge>& edges, const std::string& filename = "cell1Ds.txt");

void exportCell2Ds(const std::vector<face>& faces, const std::string& filename = "cell2Ds.txt");

void exportCell3Ds(const std::vector<polyhedron>& polyhedra, const std::string& filename = "cell3Ds.txt");

void normalize(vertex& v);

void buildTetrahedron(std::vector<vertex>& vertices, std::vector<edge>& edges, std::vector<face>& faces, polyhedron& polyhedron);

void buildEsahedron(std::vector<vertex>& vertices, std::vector<edge>& edges, std::vector<face>& faces, polyhedron& polyhedron);

void buildOctahedron(std::vector<vertex>& vertices, std::vector<edge>& edges, std::vector<face>& faces, polyhedron& polyhedron);

void buildDodecahedron(std::vector<vertex>& vertices, std::vector<edge>& edges, std::vector<face>& faces, polyhedron& polyhedron);

void buildIcosahedron(std::vector<vertex>& vertices, std::vector<edge>& edges, std::vector<face>& faces, polyhedron& polyhedron);

void buildPolyhedron(int p, int q, int b, int c, std::vector<vertex>& vertices, std::vector<edge>& edges, std::vector<face>& faces, polyhedron& polyhedron);


