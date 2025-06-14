#pragma once
#include <gtest/gtest.h>
#include "PolyhedronMesh.hpp"
#include "Utils.hpp"

using namespace std;
using namespace PolyhedronMesh;

void checkPolyhedronStructure(const vector<vertex>& vertices, const vector<edge>& edges, const vector<face>& faces, const polyhedron& poly, int expectedV, int expectedE, int expectedF) {
    EXPECT_EQ(vertices.size(), expectedV);
    EXPECT_EQ(edges.size(), expectedE);
    EXPECT_EQ(faces.size(), expectedF);

    for (const face& f : faces) {
        EXPECT_TRUE(isFaceConsistent(f, edges));
    }

    for (int vid : poly.vertex_ids) {
        EXPECT_GE(vid, 0);
        EXPECT_LT(vid, vertices.size());
    }

    for (int eid : poly.edge_ids) {
        EXPECT_GE(eid, 0);
        EXPECT_LT(eid, edges.size());
    }

    for (int fid : poly.face_ids) {
        EXPECT_GE(fid, 0);
        EXPECT_LT(fid, faces.size());
    }
}

//verifica solidi platonici
TEST(PolyhedronTest, Tetrahedron) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    buildTetrahedron(verts, edges, faces, p);
    checkPolyhedronStructure(verts, edges, faces, p, 4, 6, 4);
}

TEST(PolyhedronTest, Esahedron) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    buildEsahedron(verts, edges, faces, p);
    checkPolyhedronStructure(verts, edges, faces, p, 8, 12, 6);
}

TEST(PolyhedronTest, Octahedron) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    buildOctahedron(verts, edges, faces, p);
    checkPolyhedronStructure(verts, edges, faces, p, 6, 12, 8);
}

TEST(PolyhedronTest, Dodecahedron) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    buildDodecahedron(verts, edges, faces, p);
    checkPolyhedronStructure(verts, edges, faces, p, 20, 30, 12);
}

TEST(PolyhedronTest, Icosahedron) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    buildIcosahedron(verts, edges, faces, p);
    checkPolyhedronStructure(verts, edges, faces, p, 12, 30, 20);
}

//test sui geodesici
TEST(GeodesicPolyhedronTest, ClassI_3_4_1_0) {
    vector<vertex> vertices;
    vector<edge> edges;
    vector<face> faces;
    polyhedron poly;

    buildClassIGeodesic(3, 4, 1, vertices, edges, faces, poly);

    // Per b = 1, otteniamo direttamente l'ottaedro base
    EXPECT_EQ(vertices.size(), 6);
    EXPECT_EQ(edges.size(), 12);
    EXPECT_EQ(faces.size(), 8);
}

TEST(GeodesicPolyhedronTest, ClassI_3_4_2_0) {
    vector<vertex> vertices;
    vector<edge> edges;
    vector<face> faces;
    polyhedron poly;

    int b = 2;
    buildClassIGeodesic(3, 4, b, vertices, edges, faces, poly);
    int T = b * b; 
    EXPECT_EQ(vertices.size(), 4 * T + 2);
    EXPECT_EQ(edges.size(), 12 * T);
    EXPECT_EQ(faces.size(), 8 * T);
}



TEST(GeodesicPolyhedronTest, ClassII_3_5_1_1) {
    vector<vertex> vertices;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    
    buildClassIIGeodesic(3, 5, 1, 1, vertices, edges, faces, p);

    int N = 1 + 1;
    int expectedFaces = 20 * N * N; 
    EXPECT_EQ(faces.size(), expectedFaces);
}


TEST(GeodesicPolyhedronTest, ClassII_3_5_2_2) {
    vector<vertex> vertices;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;

    buildClassIIGeodesic(3, 5, 2, 2, vertices, edges, faces, p);

    int N = 2 + 2;
    int expectedFaces = 20 * N * N; 
    EXPECT_EQ(faces.size(), expectedFaces);
}


TEST(GeodesicPolyhedronTest, ClassII_3_3_1_1) {
    vector<vertex> vertices;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;

    buildClassIIGeodesic(3, 3, 1, 1, vertices, edges, faces, p);

    int N = 1 + 1;
    int expectedFaces = 4 * N * N; // Tetraedro ha 4 facce triangolari
    EXPECT_EQ(faces.size(), expectedFaces);
}



//test duali
TEST(DualPolyhedronTest, TetrahedronDual) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron original;
    buildTetrahedron(verts, edges, faces, original);

    vector<vertex> dualVerts;
    vector<face> dualFaces;
	vector<edge> dualEdges;
    polyhedron dual;
    buildDualPolyhedron(verts, edges, faces, original, dualVerts, dualEdges, dualFaces, dual);

    EXPECT_EQ(dualVerts.size(), faces.size());
    EXPECT_EQ(dualFaces.size(), verts.size());
}

TEST(DualPolyhedronTest, OctahedronDual) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron original;
    buildOctahedron(verts, edges, faces, original);

    vector<vertex> dualVerts;
    vector<face> dualFaces;
	vector<edge> dualEdges;
    polyhedron dual;
    buildDualPolyhedron(verts, edges, faces, original, dualVerts, dualEdges, dualFaces, dual);

    EXPECT_EQ(dualVerts.size(), faces.size());
    EXPECT_EQ(dualFaces.size(), verts.size());
}

TEST(DualPolyhedronTest, IcosahedronDual) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
	vector<edge> dualEdges;
    polyhedron original;
    buildIcosahedron(verts, edges, faces, original);

    vector<vertex> dualVerts;
    vector<face> dualFaces;
    polyhedron dual;
    buildDualPolyhedron(verts, edges, faces, original, dualVerts, dualEdges, dualFaces, dual);

    EXPECT_EQ(dualVerts.size(), faces.size());
    EXPECT_EQ(dualFaces.size(), verts.size());
}

TEST(ShortestPathTest, From0to3_Tetrahedron) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    buildTetrahedron(verts, edges, faces, p);

    findShortestPath(verts, edges, 0, 3);

    int pathVertices = 0, pathEdges = 0;
    for (auto& v : verts) pathVertices += v.ShortPath;
    for (auto& e : edges) pathEdges += e.ShortPath;

    EXPECT_GE(pathVertices, 2); 
    EXPECT_GE(pathEdges, 1);    
    EXPECT_LE(pathEdges, 2);    
}

TEST(ShortestPathTest, From0to3_Esahedron) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    buildEsahedron(verts, edges, faces, p);

    findShortestPath(verts, edges, 0, 3);

    int pathVertices = 0, pathEdges = 0;
    for (auto& v : verts) pathVertices += v.ShortPath;
    for (auto& e : edges) pathEdges += e.ShortPath;

    EXPECT_GE(pathVertices, 2);
    EXPECT_GE(pathEdges, 1);
}

TEST(ShortestPathTest, From0to4_Octahedron) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    buildOctahedron(verts, edges, faces, p);

    findShortestPath(verts, edges, 0, 4);

    int pathVertices = 0, pathEdges = 0;
    for (auto& v : verts) pathVertices += v.ShortPath;
    for (auto& e : edges) pathEdges += e.ShortPath;

    EXPECT_GE(pathVertices, 2);
    EXPECT_GE(pathEdges, 1);
}

TEST(ShortestPathTest, From0to4_Dodecahedron) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    buildDodecahedron(verts, edges, faces, p);

    findShortestPath(verts, edges, 0, 4);

    int pathVertices = 0, pathEdges = 0;
    for (auto& v : verts) pathVertices += v.ShortPath;
    for (auto& e : edges) pathEdges += e.ShortPath;

    EXPECT_GE(pathVertices, 2);
    EXPECT_GE(pathEdges, 1);
}

TEST(ShortestPathTest, From0to5_Icosahedron) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p;
    buildIcosahedron(verts, edges, faces, p);

    findShortestPath(verts, edges, 0, 5);

    int pathVertices = 0, pathEdges = 0;
    for (auto& v : verts) pathVertices += v.ShortPath;
    for (auto& e : edges) pathEdges += e.ShortPath;

    EXPECT_GE(pathVertices, 2);  
    EXPECT_GE(pathEdges, 1);     
}

