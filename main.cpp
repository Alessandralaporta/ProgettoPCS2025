#include <iostream>
#include <fstream>
#include <string>
#include "Utils.hpp"
#include "PolyhedronMesh.hpp"

int main() {
	using namespace std;
	using namespace PolyhedronMesh;
	
	vector<vertex> vertices;
    vector<edge> edges;
    vector<face> faces;
    polyhedron poly;
	//buildTetrahedron(vertices, edges, faces, poly);
    // buildEsahedron(vertices, edges, faces, poly);
    //buildOctahedron(vertices, edges, faces, poly);
    //buildDodecahedron(vertices, edges, faces, poly); 
    // buildIcosahedron(vertices, edges, faces, poly);   

	int p = 3, q = 5, b = 3, c = 0;
	if (b == 0 && c == 0) {
        buildPolyhedron(p, q, b, c, vertices, edges, faces, poly);  // diretto
    } else {
        buildClassIGeodesic(p, q, b, vertices, edges, faces, poly);  // geodesico
        //projectVerticesOnUnitSphere(vertices);
    }

    cout << "Geometria costruita: " << vertices.size() << " vertici, "
         << edges.size() << " spigoli, " << faces.size() << " facce\n";

    int startId = 0;
    int endId = vertices.size() / 2;  // metà arbitraria
    findShortestPath(vertices, edges, startId, endId);

    string outputDirectory = ".";
    exportToParaview(vertices, edges, outputDirectory);

    cout << "Esportazione completata in: " << outputDirectory << endl;
	
	/*buildDualPolyhedron(vertices, faces, original, dualVertices, dualFaces, dualPoly);
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
