//  File: topology.hpp
//  Role: Define classes required for topological analysis

#ifndef _topology_hpp_
#define _topology_hpp_

#include <vector>
#include <cassert>
#include <algorithm>
#include <utility>
#include <string>
#include <map>
#include <set>
#include <sstream>

using std::vector;
using std::pair;
using std::string;
using std::ostream;
using std::istream;
using std::set;
using std::map;
typedef unsigned short id_type;

class Grain {
public:
	// CONSTRUCTORS
	Grain(const id_type& i);
	~Grain() {}

	// ACCESSORS
	const id_type getID() const {
		return me;
	}
	float getVolume() {
		return volume;
	}
	const float getVolume() const {
		return volume;
	}
	const Point<float> getCentroid() const {
		return centroid / volume;
	}
	const int numEdges() const {
		return edges.size();
	}
	const int numFaces() const {
		return faces.size();
	}
	const int numVerts() const {
		return vertices.size();
	}
	const int getEuler() const {
		return vertices.size() + faces.size() - edges.size();
	}
	const bool isExcluded() const {
		return excluded;
	}
	template <int dim>
	const bool isPlateaunic() const {
		if (dim==2) return (vertices.size()==edges.size());
		return (edges.size() == (3 * faces.size() - 6));
	}
	template <int dim>
	const vector<int>& getPvector() const {
		return pvector;
	}
	const set<id_type>& getFaces() const {
		return faces;
	}
	const map<set<id_type>, vector<Point<int> > >& getEdges() const {
		return edges;
	}

	// MODIFIERS
	template <int dim> int numVerts();
	void addMass(float v, const Point<float>& p) {
		centroid += v * p;
		volume += v;
	}
	void addVolume(float v = 1.0) {
		volume += v;
	}
	void setExclusion(const bool b = true) {
		excluded = b;
	}
	void updateEdges(const set<id_type>& e, const Point<int>& p) {
		edges[e].push_back(p);
	}
	void updateVertices(const set<id_type>& v) {
		vertices.insert(v);
	}
	void addNeighbor(const id_type& i) {
		faces.insert(i);
	}
	template <int dim> void estimateCentroid();
	template <int dim> vector<int> computeTopology();
	template <int dim> float getPbar();

private:
	// INTERNALS
	template <int dim> int inferVertices();
	template <int dim> int checkFaces();
	template <int dim> int makePvector(vector<int>& pvec);

	//REPRESENTATION
	id_type me;
	bool excluded;
	bool face_check;
	bool vert_check;
	bool topo_check;
	Point<float> centroid;
	float volume;
	set<id_type> faces;
	vector<int> pvector;
	set<set<id_type> > vertices;
	map<set<id_type>, vector<Point<int> > > edges;

};

#endif
