#include "Utils.hpp"
#include <fstream>
#include <iostream>

using namespace std;

//scrivere gli id a partire da 0
//ordinare in cell2Ds

void exportCell0Ds (const vector<vertex>& vertices, const string& filename = "cell0Ds.txt")
{
	ofstream file(filename);
	if (!ile.is_open()){
		cerr << "Errore nell'apertura del file Cell0Ds" << filename << endl;
		return;
	}
	
	for (const auto& v:vertices){
		file << v.id << " " << v.x << " " << v.y << " " << v.z << " " << v.shortPath << endl;		
	}
	file.close();
}

void exportCell1Ds (const vector<edge>& edges, const string& filename = "cell1Ds.txt")
{
	ofstream file(filename);
	if (!ile.is_open()){
		cerr << "Errore nell'apertura del file Cell1Ds" << filename << endl;
		return;
	}
	
	for (const auto& e:edges){
		file << e.id << " " << ID.origin << " " << ID.end << endl;		
	}
	file.close();
}

void exportCell2Ds (const vector<face>& faces, const string& filename = "cell2Ds.txt")
{
	ofstream file(filename);
	if (!ile.is_open()){
		cerr << "Errore nell'apertura del file Cell2Ds" << filename << endl;
		return;
	}
	
	for (const auto& f:faces){
		file << f.id << " " << num.vertices << " " << num.edges << " " << ID.vertices << " " << ID.edges << endl;	
	}
	file.close();
}

void exportCell3Ds (const vector<polyhedron>& polyhedra, const string& filename = "cell3Ds.txt")
{
	ofstream file(filename);
	if (!ile.is_open()){
		cerr << "Errore nell'apertura del file Cell3Ds" << filename << endl;
		return;
	}
	
	for (const auto& p:polyhedra){
		file << p.id << " " << num.vertices << " " << num_edges << " " << num_faces << ID.vertices << " " << ID.edges << " " << ID.faces << endl;		
	}
	file.close();
}