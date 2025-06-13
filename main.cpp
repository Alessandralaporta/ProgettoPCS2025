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
	
	int p = 3, q = 5, b = 2, c = 0;
	
	buildPolyhedron(p, q, b, c, vertices, edges, faces, poly);

    cout << "Geometria costruita: " << vertices.size() << " vertici, "
         << edges.size() << " spigoli, " << faces.size() << " facce\n";
	
    int startId = 0;
    int endId = vertices.size() / 2;  // metà arbitraria
    findShortestPath(vertices, edges, startId, endId);

    string outputDirectory = ".";
    //exportToParaview(vertices, edges, outputDirectory);
	exportCell0Ds(vertices, "Cell0Ds.txt");
	exportCell1Ds(edges, "Cell1Ds.txt");
	exportCell2Ds(faces, "Cell2Ds.txt");
	exportCell3Ds({poly}, "Cell3Ds.txt");
	
	vector<vertex> dualVertices;
	vector<edge> dualEdges;
    vector<face> dualFaces;
    polyhedron dualPoly;

    buildDualPolyhedron(vertices, edges, faces, poly, dualVertices, dualEdges, dualFaces, dualPoly);
	cout << "Num Vertici duali: " << dualVertices.size() << std::endl;
	cout << "Num Spigoli duali: " << dualEdges.size() << std::endl;
	cout << "Num Facce duali:   " << dualFaces.size() << std::endl;
	exportToParaview(dualVertices, dualEdges, "."); 

    cout << "Esportazione completata in: " << outputDirectory << endl;
	return 0;
}
