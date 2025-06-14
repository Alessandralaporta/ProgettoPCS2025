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

    int p = 3, q = 4, b = 2, c = 2;
	//buildClassIIGeodesic(p,q,b,c,vertices,edges,faces,poly);
    buildDualFromBaseThenGeodesic(p, q, b, c, vertices, edges, faces, poly);

    int startId = 0;
    int endId = vertices.size() / 2;  // metà arbitraria
    findShortestPath(vertices, edges, startId, endId);

    exportCell0Ds(vertices, "Cell0Ds.txt");
    exportCell1Ds(edges, "Cell1Ds.txt");
    exportCell2Ds(faces, "Cell2Ds.txt");
    exportCell3Ds({poly}, "Cell3Ds.txt");

    exportToParaview(vertices, edges, ".");

    return 0;
}
