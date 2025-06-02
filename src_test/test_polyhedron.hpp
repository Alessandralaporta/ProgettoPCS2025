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
void testGeodesicCounts(int p, int q, int b, int c) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron p_out;
    buildGeodesicPolyhedron(p, q, b, c, verts, edges, faces, p_out);

    int T = b*b + b*c + c*c;
    if (q == 3) {
        EXPECT_EQ(verts.size(), 2*T + 2);
        EXPECT_EQ(edges.size(), 6*T);
        EXPECT_EQ(faces.size(), 4*T);
    } else if (q == 4) {
        EXPECT_EQ(verts.size(), 4*T + 2);
        EXPECT_EQ(edges.size(), 12*T);
        EXPECT_EQ(faces.size(), 8*T);
    } else if (q == 5) {
        EXPECT_EQ(verts.size(), 10*T + 2);
        EXPECT_EQ(edges.size(), 30*T);
        EXPECT_EQ(faces.size(), 20*T);
    }
}

//dovrebbe funzionare anche con ClassI_3_5_3_0 e ClassII_3_5_3_3 e altri valori, da verificare!
//e ovviamente dovrebbe funzionare anche per tutti i solidi
TEST(GeodesicPolyhedronTest, ClassI_3_5_1_0) {
    testGeodesicCounts(3, 5, 1, 0);
}

TEST(GeodesicPolyhedronTest, ClassII_3_5_1_1) {
    testGeodesicCounts(3, 5, 1, 1);
}

TEST(GeodesicPolyhedronTest, ClassI_3_5_2_0) {
    testGeodesicCounts(3, 5, 2, 0);
}

TEST(GeodesicPolyhedronTest, ClassII_3_5_2_2) {
    testGeodesicCounts(3, 5, 2, 2);
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
    polyhedron dual;
    buildDualPolyhedron(verts, faces, original, dualVerts, dualFaces, dual);

    EXPECT_EQ(dualVerts.size(), faces.size());
    EXPECT_EQ(dualFaces.size(), verts.size());
}

TEST(DualPolyhedronTest, EsahedronDual) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron original;
    buildEsahedron(verts, edges, faces, original);

    vector<vertex> dualVerts;
    vector<face> dualFaces;
    polyhedron dual;
    buildDualPolyhedron(verts, faces, original, dualVerts, dualFaces, dual);

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
    polyhedron dual;
    buildDualPolyhedron(verts, faces, original, dualVerts, dualFaces, dual);

    EXPECT_EQ(dualVerts.size(), faces.size());
    EXPECT_EQ(dualFaces.size(), verts.size());
}

TEST(DualPolyhedronTest, DodecahedronDual) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron original;
    buildDodecahedron(verts, edges, faces, original);

    vector<vertex> dualVerts;
    vector<face> dualFaces;
    polyhedron dual;
    buildDualPolyhedron(verts, faces, original, dualVerts, dualFaces, dual);

    EXPECT_EQ(dualVerts.size(), faces.size());
    EXPECT_EQ(dualFaces.size(), verts.size());
    EXPECT_EQ(dual.vertex_ids.size(), dualVerts.size());
    EXPECT_EQ(dual.face_ids.size(), dualFaces.size());
}

TEST(DualPolyhedronTest, IcosahedronDual) {
    vector<vertex> verts;
    vector<edge> edges;
    vector<face> faces;
    polyhedron original;
    buildIcosahedron(verts, edges, faces, original);

    vector<vertex> dualVerts;
    vector<face> dualFaces;
    polyhedron dual;
    buildDualPolyhedron(verts, faces, original, dualVerts, dualFaces, dual);

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

//GeodesicPolyhedronTest