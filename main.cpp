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

    int p = 3, q = 5, b = 3, c = 3;

    // A questo punto possiamo usare la funzione centralizzata che si occupa di calcolare Class I e Class II
    buildClassIIGeodesic(p, q, b, c, vertices, edges, faces, poly);

    // Ricerca del cammino più breve (arbitrario, ad esempio tra il primo e l'ultimo vertice)
    int startId = 0;
    int endId = vertices.size() / 2;  // metà arbitraria
    findShortestPath(vertices, edges, startId, endId);

    // Esportazione dei dati
    exportCell0Ds(vertices, "Cell0Ds.txt");
    exportCell1Ds(edges, "Cell1Ds.txt");
    exportCell2Ds(faces, "Cell2Ds.txt");
    exportCell3Ds({poly}, "Cell3Ds.txt");

    exportToParaview(vertices, edges, ".");

    cout << "Esportazione completata in: " << "." << endl;

    return 0;
}
