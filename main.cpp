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
	
	int p = 5, q = 3, b = 0, c = 0;
	
    bool classI = (b > 0 && c == 0) || (b == 0 && c > 0);
    bool classII = (b > 0 && b == c);
    if (b == 0 && c == 0) {

        buildPolyhedron(p, q, b, c, vertices, edges, faces, poly);
        cout << "Costruzione poliedro base...\n";
    } else if (classI) {
        cout << "Costruzione geodesico Classe I (b = " << b << ", c = " << c << ")...\n";
        buildClassIGeodesic(p, q, max(b, c), vertices, edges, faces, poly);
    } else if (classII) {
        cout << "Costruzione geodesico Classe II (b = c = " << b << ")...\n";
        buildClassIIGeodesic(p, q, b, c, vertices, edges, faces, poly);
    } else {
        cerr << "Parametri non validi o classe non supportata!\n";
        return 1;
    }

    cout << "Geometria costruita: " << vertices.size() << " vertici, "
         << edges.size() << " spigoli, " << faces.size() << " facce\n";

    int startId = 0;
    int endId = vertices.size() / 2;  // metà arbitraria
    findShortestPath(vertices, edges, startId, endId);

    string outputDirectory = ".";
    exportToParaview(vertices, edges, outputDirectory);
	exportCell0Ds(vertices, "Cell0Ds.txt");
	exportCell1Ds(edges, "Cell1Ds.txt");
	exportCell2Ds(faces, "Cell2Ds.txt");
	exportCell3Ds({poly}, "Cell3Ds.txt");
	

    cout << "Esportazione completata in: " << outputDirectory << endl;
	
	/*buildDualPolyhedron(vertices, faces, original, dualVertices, dualFaces, dualPoly);
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
