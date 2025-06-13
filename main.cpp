#include <iostream>
#include <fstream>
#include <string>
#include "Utils.hpp"
#include "PolyhedronMesh.hpp"
#include <set>
int main() {
	using namespace std;
	using namespace PolyhedronMesh;
	
	vector<vertex> vertices;
    vector<edge> edges;
    vector<face> faces;
    polyhedron poly;
	
	int p = 3, q = 3, b = 3, c = 0;
	
	buildDualFromBaseThenGeodesic(p, q, b, c, vertices, edges, faces, poly);

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

	cout << "Num Vertici duali: " << vertices.size() << std::endl;
	cout << "Num Spigoli duali: " << edges.size() << std::endl;
	cout << "Num Facce duali:   " << faces.size() << std::endl;
	exportToParaview(vertices, edges, "."); 

    cout << "Esportazione completata in: " << outputDirectory << endl;
	return 0;
}
