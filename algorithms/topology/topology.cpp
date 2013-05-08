//  File: topology.cpp
//  Role: Implement functions required to analyze grain topology.

#include <iostream>
#include <algorithm>
#include "topology.hpp"

Grain::Grain(const id_type& i) {
	me = i;
	centroid = Point<float>(0., 0., 0.);
	excluded = false;   // is the grain excluded from analysis?
	vert_check = false; // have vertices been inferred yet?
	face_check = false; // have faces been checked yet?
	volume = 0.0;
}

void Grain::estimateCentroid() {
	if (faces.size() == 0 || edges.size() == 0) return;
	int count = 0;
	Point<int> min(1000, 1000, 1000);
	Point<int> max(0, 0, 0);
	// The grain has no concept of the domain size. Sweep through once to bound the domain.
	for (map<set<id_type>, vector<Point<int> > >::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
		count += itr->second.size();
		for (int i = 0; i < itr->second.size(); ++i) {
			Point<int> x = itr->second[i];
			if (x.x < min.x) min.x = x.x;
			else if (x.x > max.x) max.x = x.x;
			else if (x.x < min.x) min.x = x.x;
			else if (x.y > max.y) max.y = x.y;
			else if (x.y < min.y) min.y = x.y;
			else if (x.z < min.z) min.z = x.z;
			else if (x.z > max.z) max.z = x.z;
		}
	}
	float weight = volume / float(count); // Final volume should match initial volume.
	// Reset centroid. Estimate it.
	volume = 0.0f;
	centroid = Point<float>(0., 0., 0.);
	for (map<set<id_type>, vector<Point<int> > >::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
		for (int i = 0; i < itr->second.size(); ++i) {
			Point<int> x = itr->second[i];
			if (max.x - x.x < x.x - min.x) x.x = x.x - max.x; //point is closer to upper than lower boundary
			if (max.y - x.y < x.y - min.y) x.y = x.y - max.y; //point is closer to upper than lower boundary
			if (max.z - x.z < x.z - min.z) x.z = x.z - max.z; //point is closer to upper than lower boundary
			this->addMass(weight, x);
		}
	}
}

vector<int> Grain::computeTopology() {
	vector<int> answer;
	answer.push_back(this->checkFaces());
	answer.push_back(this->inferVertices());
	answer.push_back(this->makePvector(pvector));
	answer.push_back(this->numFaces());
	return answer;
}

int Grain::numVerts() {
	if (edges.size() == 0) return 0;
	if (!vert_check) this->inferVertices();
	return vertices.size();
}

int Grain::inferVertices() {
	vert_check = true;
	vertices.clear();
	if (this->numEdges() < 3) return -1;
	if (this->numFaces() < 3) return 0;
	for (map<set<id_type>, vector<Point<int> > >::iterator eitr = edges.begin(); eitr != edges.end(); ++eitr) {
		id_type a = *(eitr->first.begin());
		id_type b = *(++eitr->first.begin());
		assert(a != b);
		set<id_type> pair_a;
		set<id_type> pair_b;
		for (set<id_type>::iterator fitr = faces.begin(); fitr != faces.end(); ++fitr) {
			id_type c = *fitr;
			if ((c == a) || (c == b)) continue;
			pair_a.clear();
			pair_a.insert(a);
			pair_a.insert(c);
			pair_b.clear();
			pair_b.insert(b);
			pair_b.insert(c);
			if ((edges.find(pair_a) != edges.end()) && (edges.find(pair_b) != edges.end())) {
				id_type grains[] = {a, b, c};
				std::sort(grains, grains + 3);
				vertices.insert(set<id_type>(grains, grains + 3));
			}
		}
	}
	return vertices.size();
}

int Grain::checkFaces() {
	//if (!vert_check) this->inferVertices();
	// Make sure that each neighbor shows up more than once in the map of edges
	int answer = 0;
	while (face_check == false) {
		face_check = true;
		for (set<id_type>::iterator fitr = faces.begin(); fitr != faces.end(); ++fitr) {
			int p = 0;
			for (map<set<id_type>, vector<Point<int> > >::iterator eitr = edges.begin(); eitr != edges.end(); ++eitr) {
				for (set<id_type>::iterator set_itr = eitr->first.begin(); set_itr != eitr->first.end(); ++set_itr) {
					if (*set_itr == *fitr) ++p; // add an edge to the specified face
				}
			}
			// Remove faces with one edge!
			if (p == 1) {
				face_check = false;
				vert_check = false;
				++answer;
				for (map<set<id_type>, vector<Point<int> > >::iterator eitr = edges.begin(); eitr != edges.end(); ++eitr) {
					for (set<id_type>::iterator set_itr = eitr->first.begin(); set_itr != eitr->first.end(); ++set_itr) {
						// Remove edge with offending neighbor
						if ((*set_itr == *fitr) && (eitr != edges.end())) {
							edges.erase( eitr );
							eitr = edges.begin();
							if (edges.size() == 0) break;
							set_itr = eitr->first.begin();
							break; // restart outer loop: check for more invalid faces
						}
						if (edges.size() == 0) break;
					}
					if (edges.size() == 0) break;
				}
				set<id_type>::iterator temp_itr(fitr);
				--fitr;
				faces.erase(temp_itr);
			}
		}
	}
	return answer;
}

int Grain::makePvector(vector<int>& pvec) {
	if (!face_check) this->checkFaces();
	if (!vert_check) this->inferVertices();
	pvec.clear();
	pvec.push_back( 0 );
	if (this->numEdges() == 0) return 0;
	for (set<id_type>::iterator fitr = faces.begin(); fitr != faces.end(); ++fitr) {
		int p = 0;
		for (map<set<id_type>, vector<Point<int> > >::iterator eitr = edges.begin(); eitr != edges.end(); ++eitr) {
			for (set<id_type>::iterator set_itr = eitr->first.begin(); set_itr != eitr->first.end(); ++set_itr) {
				if (*set_itr == *fitr) ++p;
			}
		}
		while (pvec.size() <= p) pvec.push_back(0);
		++pvec[p];
	}
	assert(pvec[1] == 0);
	assert(pvec.size() < 2 * (this->numFaces()));
	return pvec.size();
}

float Grain::getPbar() {
	if (!face_check || !vert_check) {
		this->checkFaces();
	}
	if (pvector.size() == 0) this->makePvector(pvector);
	int answer = 0, count = 0;
	for (int i = 0; i < pvector.size(); ++i) {
		count += pvector[i]; // number of faces
		answer += pvector[i] * i; // number of edges
	}
	return float(answer) / count;
}

void printCSV(std::ofstream& o, const Grain& g, const int& p, const int& N) {
	if ( g.isExcluded() && g.numFaces() > 0 ) {
		Point<int> x = g.getCentroid();
		o << g.getID() << ',' << x.x << ',' << x.y << ',' << x.z << ',' << g.getVolume() << ",0," << g.numFaces() << ",0,0";
		const vector<int>& pvec = g.getPvector();
		assert(pvec.size() != 0);
		for (int i = 2; i < p; ++i) o << ",0";
		const set<id_type>& faces = g.getFaces();
		for (set<id_type>::iterator itr = faces.begin(); itr != faces.end(); ++itr) o << ',' << *itr;
		for (int i = faces.size(); i < N; ++i) o << ',';
		o << '\n';
	} else if ( g.numVerts() < 2 || g.numEdges() < 3 || g.numFaces() < 3 ) return; // Fewer F, V, or E than a Brazil nut? Invalid.
	else {
		Point<int> x = g.getCentroid();
		o << g.getID() << ',' << x.x << ',' << x.y << ',' << x.z << ',' << g.getVolume() << ',' << g.numVerts() << ',' << g.numFaces() << ',' << g.numEdges() << ',' << g.getEuler();
		const vector<int>& pvec = g.getPvector();
		assert(pvec.size() != 0);
		for (int i = 2; i < pvec.size(); ++i) o << ',' << pvec[i];
		for (int i = pvec.size(); i < p; ++i) o << ",0";
		const set<id_type>& faces = g.getFaces();
		for (set<id_type>::iterator itr = faces.begin(); itr != faces.end(); ++itr) o << ',' << *itr;
		for (int i = faces.size(); i < N; ++i) o << ',';
		o << '\n';
	}
}
