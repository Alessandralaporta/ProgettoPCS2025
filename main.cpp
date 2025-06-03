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
	//buildTetrahedron(vertices, edges, faces, poly);
    // buildEsahedron(vertices, edges, faces, poly);
    buildOctahedron(vertices, edges, faces, poly);
    // buildDodecahedron(vertices, edges, faces, poly); // ⚠️ ancora da correggere
    // buildIcosahedron(vertices, edges, faces, poly);   // ⚠️ ancora da correggere

    cout << "Verifica consistenza delle facce del poliedro: " << poly.id << "\n";

    for (const auto& f : faces) {
        bool consistent = isFaceConsistent(f, edges, 1e-6);
        cout << "Faccia " << f.id << ": " << (consistent ? "✅ OK" : "❌ INCONSISTENTE") << "\n";
    }
	/*
	int p = 3, q = 5, b = 2, c = 2; // esempio: classe I cubo
	if (b == 0 && c == 0) {
        buildPolyhedron(p, q, b, c, vertices, edges, faces, poly);  // diretto
    } else {
        buildGeodesicPolyhedron(p, q, b, c, vertices, edges, faces, poly);  // geodesico
        projectVerticesOnUnitSphere(vertices);
    }
	
	//buildDualPolyhedron(vertices, faces, original, dualVertices, dualFaces, dualPoly);
	exportCell0Ds(vertices, "Cell0Ds.txt");
	exportCell1Ds(edges, "Cell1Ds.txt");
	exportCell2Ds(faces, "Cell2Ds.txt");
	exportCell3Ds({poly}, "Cell3Ds.txt");
	
	vector<vertex> vertices;
    vector<edge> edges;
    vector<face> faces;
    polyhedron poly;

    buildEsahedron(vertices, edges, faces, poly);

    vector<vertex> dualVertices;
    vector<face> dualFaces;
    polyhedron dualPoly;

    buildDualPolyhedron(vertices, faces, poly, dualVertices, dualFaces, dualPoly);

    vector<polyhedron> polyhedra = {poly, dualPoly};
    exportCell3Ds(polyhedra, "Cell3Ds.txt");*/
	return 0;
}
