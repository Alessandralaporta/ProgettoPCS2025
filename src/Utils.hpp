#pragma once
#include <string>
#include <vector>
#include "PolyhedronMesh.hpp"
#include "UCDUtilities.hpp"

void normalize(PolyhedronMesh::vertex& v);

void buildTetrahedron(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

void buildEsahedron(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

void buildOctahedron(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

void buildDodecahedron(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

void buildIcosahedron(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

void buildPolyhedron(int p, int q, int b, int c, std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& polyhedron);

bool sameVertex(const PolyhedronMesh::vertex& a, const PolyhedronMesh::vertex& b, double tolerance = 1e-6);

bool isFaceConsistent(const PolyhedronMesh::face& face, const std::vector<PolyhedronMesh::edge>& edges, double tolerance = 1e-4);

int getOrAddVertex(double x, double y, double z, std::vector<PolyhedronMesh::vertex>& vertices);

int getOrAddVertex(double x, double y, double z, std::vector<PolyhedronMesh::vertex>& vertices, double tolerance);

void triangulateFaces(const std::vector<PolyhedronMesh::face>& inputFaces, std::vector<PolyhedronMesh::face>& outputFaces);

void buildClassIGeodesic(int p, int q, int b, std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& poly);

void buildClassIIGeodesic(int p, int q, int b, int c, std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, std::vector<PolyhedronMesh::face>& faces, PolyhedronMesh::polyhedron& poly);

double distance(const PolyhedronMesh::vertex& a, const PolyhedronMesh::vertex& b);

void findShortestPath(std::vector<PolyhedronMesh::vertex>& vertices, std::vector<PolyhedronMesh::edge>& edges, int startId, int endId);

void exportCell0Ds(const std::vector<PolyhedronMesh::vertex>& vertices, const std::string& filename = "cell0Ds.txt");

void exportCell1Ds(const std::vector<PolyhedronMesh::edge>& edges, const std::string& filename = "cell1Ds.txt");

void exportCell2Ds(const std::vector<PolyhedronMesh::face>& faces, const std::string& filename = "cell2Ds.txt");

void exportCell3Ds(const std::vector<PolyhedronMesh::polyhedron>& polyhedra, const std::string& filename = "cell3Ds.txt");

void projectVerticesOnUnitSphere(std::vector<PolyhedronMesh::vertex> & vertices);

std::vector<PolyhedronMesh::vertex> calculateCentroids(const std::vector<PolyhedronMesh::vertex>& vertices, const std::vector<PolyhedronMesh::face>& faces);

void buildDualPolyhedron(const std::vector<PolyhedronMesh::vertex>& vertices, const std::vector<PolyhedronMesh::face>& faces, const PolyhedronMesh::polyhedron& original, std::vector<PolyhedronMesh::vertex>& dualVertices, std::vector<PolyhedronMesh::face>& dualFaces, PolyhedronMesh::polyhedron& dualPoly);

void exportToParaview(const std::vector<PolyhedronMesh::vertex>& vertices, const std::vector<PolyhedronMesh::edge>& edges, const std::string& outputDirectory);