#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <gtest/gtest.h>
#include "PolyhedronMesh.hpp"
#include "Utils.hpp"

using namespace std;
using namespace PolyhedronMesh;

//test caso generale normalized
TEST(UtilsTest, NormalizeRegularVector) {
    vertex v{0, 3.0, 0.0, 4.0};
    normalize(v);
    EXPECT_NEAR(v.x, 0.6, 1e-6);
    EXPECT_NEAR(v.y, 0.0, 1e-6);
    EXPECT_NEAR(v.z, 0.8, 1e-6);
}
// test caso vettore nullo normalized
TEST(UtilsTest, NormalizeZeroVector) {
    vertex v{0, 0.0, 0.0, 0.0};
    normalize(v);
    EXPECT_DOUBLE_EQ(v.x, 0.0);
    EXPECT_DOUBLE_EQ(v.y, 0.0);
    EXPECT_DOUBLE_EQ(v.z, 0.0);
}
// test vettore già normalizzato
TEST(UtilsTest, NormalizeAlreadyNormalized) {
    vertex v{0, 0.0, 0.6, 0.8};
    normalize(v);
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);;
    EXPECT_NEAR(len, 1.0, 1e-6);
    EXPECT_NEAR(v.y, 0.6, 1e-6);
    EXPECT_NEAR(v.z, 0.8, 1e-6);
}
//test componenti negative normalizzate
TEST(UtilsTest, NormalizeNegativeComponents) {
    vertex v{0, -3.0, 0.0, -4.0};
    normalize(v);
    EXPECT_NEAR(v.x, -0.6, 1e-6);
    EXPECT_NEAR(v.y, 0.0, 1e-6);
    EXPECT_NEAR(v.z, -0.8, 1e-6);
}

// test caso generale distance
TEST(UtilsTest, DistanceSimple) {
    vertex a{0, 0.0, 0.0, 0.0};
    vertex b{1, 1.0, 0.0, 0.0};
    EXPECT_NEAR(distance(a, b), 1.0, 1e-6);
}

// test distanza con lo stesso punto (distanza nulla)
TEST(UtilsTest, DistanceZero){
    vertex a{0, 1.0, 2.0, 3.0};
    EXPECT_NEAR(distance(a, a), 0.0, 1e-6);
}

//test distanza coordinate negative
TEST(UtilsTest, DistanceNegativeCoordinates) {
    vertex a{0, -1.0, -2.0, -3.0};
    vertex b{1, 4.0, 6.0, 8.0};
    double expected = sqrt((a.x - b.x) * (a.x - b.x) +  (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
    EXPECT_NEAR(distance(a, b), expected, 1e-6);
}

//test distanaza vertici identici
TEST(UtilsTest, DistanceIdenticalVertices) {
    vertex a{0, 1.0, 2.0, 3.0};
    vertex b{1, 1.0, 2.0, 3.0};
    EXPECT_DOUBLE_EQ(distance(a, b), 0.0);
}

//verifica aggiunta primo vertice inserito corrrettamente (GetOrAddVertex)
TEST(UtilsTest, GetOrAddVertexAdd) {
    vector<vertex> verts;
    int id = getOrAddVertex(1.00001, 2.00001, 3.00001, verts);
    EXPECT_EQ(id, 0);
    EXPECT_EQ(verts.size(), 1);
}

//test dopo aver aggiunto un vertice prova ad aggiungerne un altro molto vicino
TEST(UtilsTest, GetOrAddVertexReuseDefaultTolerance) {
    vector<vertex> verts;
    getOrAddVertex(1.0, 2.0, 3.0, verts);
    int id = getOrAddVertex(1.00001, 2.00001, 3.00001, verts);
    EXPECT_EQ(id, 0);
    EXPECT_EQ(verts.size(), 1);
}

// test aggiunta nuovo vertice quando è fuori dalla tolleranza
TEST(UtilsTest, GetOrAddVertexNewWhenOutsideTolerance) {
    vector<vertex> verts;
    getOrAddVertex(1.0, 2.0, 3.0, verts);
    int id = getOrAddVertex(1.01, 2.0, 3.0, verts);
    EXPECT_EQ(id, 1);
    EXPECT_EQ(verts.size(), 2);
}

// test aggiunta multipla di vertici distinti con ID corretti e dimensione giusta
TEST(UtilsTest, GetOrAddVertexMultipleAdds) {
    vector<vertex> verts;
    int id0 = getOrAddVertex(0.0, 0.0, 0.0, verts);
    int id1 = getOrAddVertex(1.0, 0.0, 0.0, verts);
    int id2 = getOrAddVertex(0.0, 1.0, 0.0, verts);
    EXPECT_EQ(id0, 0);
    EXPECT_EQ(id1, 1);
    EXPECT_EQ(id2, 2);
    EXPECT_EQ(verts.size(), 3);
}

//test costruzione vertici
TEST(UtilsTest, ProjectVerticesBasic) {
    vector<vertex> verts = {
        {0, 3.0, 0.0, 4.0},
        {1, 0.0, 0.0, 0.0}
    };
    projectVerticesOnUnitSphere(verts);
    EXPECT_NEAR(std::sqrt(verts[0].x * verts[0].x + verts[0].y * verts[0].y + verts[0].z * verts[0].z), 1.0, 1e-6); 
    EXPECT_DOUBLE_EQ(verts[1].x, 0.0);
	EXPECT_DOUBLE_EQ(verts[1].y, 0.0);
    EXPECT_DOUBLE_EQ(verts[1].z, 0.0);
}

//test costruzione vertici già normalizzati
TEST(UtilsTest, ProjectVerticesAlreadyNormalized) {
    vector<vertex> verts = {
        {0, 0.0, 0.6, 0.8}
    };
    projectVerticesOnUnitSphere(verts);
    double len = sqrt(verts[0].x*verts[0].x + verts[0].y*verts[0].y + verts[0].z*verts[0].z);
    EXPECT_NEAR(len, 1.0, 1e-6);
    EXPECT_NEAR(verts[0].y, 0.6, 1e-6);
    EXPECT_NEAR(verts[0].z, 0.8, 1e-6);
}

// test costruzione vertici con componenti negative
TEST(UtilsTest, ProjectVerticesNegativeComponents) {
    vector<vertex> verts = {
        {0, -3.0, 0.0, -4.0}
    };

    projectVerticesOnUnitSphere(verts);
    double len = sqrt(verts[0].x*verts[0].x + verts[0].y*verts[0].y + verts[0].z*verts[0].z);
    EXPECT_NEAR(len, 1.0, 1e-6);
    EXPECT_NEAR(verts[0].x, -0.6, 1e-6);
    EXPECT_NEAR(verts[0].z, -0.8, 1e-6);
}

 // test costruzione vertici con componenti molto piccole(quasi zero)
TEST(UtilsTest, ProjectVerticesVerySmallVector) {
    vector<vertex> verts = {
        {0, 1e-10, 0.0, 1e-10}
    };
    projectVerticesOnUnitSphere(verts);
    EXPECT_NEAR(verts[0].x, 0.0, 1e-8);
    EXPECT_NEAR(verts[0].y, 0.0, 1e-8);
    EXPECT_NEAR(verts[0].z, 0.0, 1e-8);
}

 // test costruzione vertici con componenti molto grandi
TEST(UtilsTest, ProjectVerticesLargeVector) {
    vector<vertex> verts = {
        {0, 1e10, 0.0, 1e10}
    };
    projectVerticesOnUnitSphere(verts);
    double len = sqrt(verts[0].x*verts[0].x + verts[0].y*verts[0].y + verts[0].z*verts[0].z);
    EXPECT_NEAR(len, 1.0, 1e-6);
}

//test samevertex
TEST(UtilsTest, ExactMatch) {
    vertex a{0, 1.0, 2.0, 3.0};
    vertex b{1, 1.0, 2.0, 3.0};
    EXPECT_TRUE(sameVertex(a, b, 0.0));
}

TEST(UtilsTest, WithinTolerance) {
    vertex a{0, 1.0, 2.0, 3.0};
    vertex b{1, 1.0005, 2.0005, 3.0005};
    EXPECT_TRUE(sameVertex(a, b, 0.001));
}

TEST(UtilsTest, OutsideToleranceX) {
    vertex a{0, 1.0, 2.0, 3.0};
    vertex b{1, 1.002, 2.0, 3.0};
    EXPECT_FALSE(sameVertex(a, b, 0.001));
}

TEST(UtilsTest, OutsideToleranceY) {
    vertex a{0, 1.0, 2.0, 3.0};
    vertex b{1, 1.0, 2.002, 3.0};
    EXPECT_FALSE(sameVertex(a, b, 0.001));
}

TEST(UtilsTest, OutsideToleranceZ) {
    vertex a{0, 1.0, 2.0, 3.0};
    vertex b{1, 1.0, 2.0, 3.002};
    EXPECT_FALSE(sameVertex(a, b, 0.001));
}
// vedere se è da mettere qui o in polyhedronMesh
// test se gli spigoli sono coerenti
// test facce meno di 3 spigoli
// test facce meno di 3 spigoli
// test facce meno di 3 spigoli
TEST(UtilsTest, LessThan3Edges) {
    face f;
    f.edge_ids = {0, 1}; // solo 2 edges

    std::vector<edge> edges;
    edges.push_back(edge(0, 0, 1)); // id=0, origin=0, end=1
    edges.push_back(edge(1, 1, 2)); // id=1, origin=1, end=2

    EXPECT_FALSE(isFaceConsistent(f, edges));
}

// test ID spigolo fuori intervallo
TEST(UtilsTest, EdgeIdOutOfRange) {
    face f;
    f.edge_ids = {0, 2, 3};

    std::vector<edge> edges;
    edges.push_back(edge(0, 0, 1));
    edges.push_back(edge(1, 1, 2));

    EXPECT_FALSE(isFaceConsistent(f, edges));
}

// test spigoli non collegati tra loro
TEST(UtilsTest, NonConsecutiveEdges) {
    face f;
    f.edge_ids = {0, 1, 2};

    std::vector<edge> edges;
    edges.push_back(edge(0, 0, 1));
    edges.push_back(edge(1, 2, 3));
    edges.push_back(edge(2, 4, 5));

    EXPECT_FALSE(isFaceConsistent(f, edges));
}

// test triangolo coerente
TEST(UtilsTest, ConsistentTriangle) {
    face f;
    f.edge_ids = {0, 1, 2};

    std::vector<edge> edges;
    edges.push_back(edge(0, 0, 1));
    edges.push_back(edge(1, 1, 2));
    edges.push_back(edge(2, 2, 0));

    EXPECT_TRUE(isFaceConsistent(f, edges));
}

// test quadrato coerente
TEST(FaceConsistencyTest, ConsistentSquare) {
    face f;
    f.edge_ids = {0, 1, 2, 3};

    std::vector<edge> edges;
    edges.push_back(edge(0, 0, 1));
    edges.push_back(edge(1, 1, 2));
    edges.push_back(edge(2, 2, 3));
    edges.push_back(edge(3, 3, 0));

    EXPECT_TRUE(isFaceConsistent(f, edges));
}

// test spigoli invertiti
TEST(UtilsTest, InvertedEdgeBreaksConsistency) {
    face f;
    f.edge_ids = {0, 1, 2};

    std::vector<edge> edges;
    edges.push_back(edge(0, 0, 1));
    edges.push_back(edge(1, 2, 1)); // invertito
    edges.push_back(edge(2, 2, 0));

    EXPECT_FALSE(isFaceConsistent(f, edges));
}
