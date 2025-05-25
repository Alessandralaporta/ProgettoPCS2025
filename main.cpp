#include <iostream>
#include <fstream>
#include <string>
#include "Utils.hpp"
#include "PolyhedronMesh.hpp"

using namespace std;
using namespace PolyhedronMesh;

int main() {
	vector<vertex> vertices;
	vector<edge> edges;
	vector<face> faces;
	polyhedron poly;

	int p = 3, q = 5, b = 2, c = 0; // esempio: classe I icosaedrica
	buildGeodesicPolyhedron(p, q, b, c, vertices, edges, faces, poly);
	exportCell0Ds(vertices, "Cell0Ds.txt");
	exportCell1Ds(edges, "Cell1Ds.txt");
	exportCell2Ds(faces, "Cell2Ds.txt");
	exportCell3Ds({poly}, "Cell3Ds.txt");
	return 0;
}