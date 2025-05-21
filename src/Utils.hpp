#pragma once
#include <string>
#include <vector>
#include "PolyhedronMesh.hpp"

void exportCell0Ds(const std::vector<PolyhedronMesh::vertex>& vertices, const std::string& filename = "cell0Ds.txt");

void exportCell1Ds(const std::vector<PolyhedronMesh::edge>& edges, const std::string& filename = "cell1Ds.txt");

void exportCell2Ds(const std::vector<PolyhedronMesh::face>& faces, const std::string& filename = "cell2Ds.txt");

void exportCell3Ds(const std::vector<PolyhedronMesh::polyhedron>& polyhedra, const std::string& filename = "cell3Ds.txt");

void normalize(PolyhedronMesh::vertex& v);

void buildTetrahedron(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

void buildEsahedron(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

void buildOctahedron(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

void buildDodecahedron(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

void buildIcosahedron(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

void buildPolyhedron(int p, int q, int b, int c, std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);


