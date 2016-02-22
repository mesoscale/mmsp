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

template <int dim>
void Grain::estimateCentroid() {
	if (edges.size() == 0) return;
	if (dim==3 && faces.size() == 0) return;
	int count = 0;
	Point<int> min(1000, 1000);
	if (dim==3) min.z=1000;
	Point<int> max(0, 0, 0);
	// The grain has no concept of the domain size. Sweep through once to bound the domain.
	for (map<set<id_type>, vector<Point<int> > >::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
		count += itr->second.size();
		for (unsigned int i = 0; i < itr->second.size(); ++i) {
			Point<int> x = itr->second[i];
			for (int d=0; d<dim; ++d) {
				if (x[d]<min[d]) min[d]=x[d];
				if (x[d]>max[d]) max[d]=x[d];
			}
		}
	}
	float weight = volume / float(count); // Final volume should match initial volume.
	// Reset centroid. Estimate it.
	volume = 0.0f;
	centroid = Point<float>(0., 0., 0.);
	for (map<set<id_type>, vector<Point<int> > >::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
		for (unsigned int i = 0; i < itr->second.size(); ++i) {
			Point<int> x = itr->second[i];
			for (int d=0; d<dim; ++d)
				if (max[d] - x[d] < x[d] - min[d]) x[d] = x[d] - max[d]; //point is closer to upper than lower boundary
			this->addMass(weight, x);
		}
	}
}

template <int dim>
vector<int> Grain::computeTopology() {
	vector<int> answer;
	if (dim==2) {
		answer.push_back(this->checkFaces<2>());
		answer.push_back(this->inferVertices<2>());
		answer.push_back(this->makePvector<2>(pvector));
		answer.push_back(this->numFaces());
	} else {
		answer.push_back(this->checkFaces<3>());
		answer.push_back(this->inferVertices<3>());
		answer.push_back(this->makePvector<3>(pvector));
		answer.push_back(this->numFaces());
	}
	return answer;
}

template <int dim>
int Grain::numVerts() {
	if (edges.size() == 0) return 0;
	if (dim==2) {
		if (!vert_check) this->inferVertices<2>();
	} else {
		if (!vert_check) this->inferVertices<3>();
	}
	return vertices.size();
}

template <int dim>
int Grain::inferVertices() {
	vert_check = true;
	if (dim==2) {
		if (this->numEdges() < 2) return -1;
		return vertices.size();
	} else {
		if (this->numEdges() < 3) return -1;
		if (this->numFaces() < 3) return 0;
		vertices.clear();
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
	}
	return vertices.size();
}

template <int dim>
int Grain::checkFaces() {
	int answer = 0;
	face_check = true;
	//if (!vert_check) this->inferVertices();
	// Make sure that each neighbor shows up more than once in the map of edges
	if (dim==2) {
		return answer;
	}
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

template <int dim>
int Grain::makePvector(vector<int>& pvec) {
	pvec.clear();
	pvec.push_back( 0 );
	if (dim==2) {
		pvec.push_back( 0 );
		if (!face_check) this->checkFaces<2>();
		if (!vert_check) this->inferVertices<2>();
		pvec.push_back(edges.size());
		return pvec.size();
	}
	if (!face_check) this->checkFaces<3>();
	if (!vert_check) this->inferVertices<3>();
	if (this->numEdges() == 0) return 0;
	for (set<id_type>::iterator fitr = faces.begin(); fitr != faces.end(); ++fitr) {
		unsigned int p = 0;
		for (map<set<id_type>, vector<Point<int> > >::iterator eitr = edges.begin(); eitr != edges.end(); ++eitr) {
			for (set<id_type>::iterator set_itr = eitr->first.begin(); set_itr != eitr->first.end(); ++set_itr) {
				if (*set_itr == *fitr) ++p;
			}
		}
		while (pvec.size() <= p) pvec.push_back(0);
		++pvec[p];
	}
	assert(pvec[1] == 0);
	return pvec.size();
}

template <int dim>
float Grain::getPbar() {
	if (dim==2) {
		if (!face_check || !vert_check)	this->checkFaces<2>();
		if (pvector.size() == 0) this->makePvector<2>(pvector);
		return 0.0f;
	}
	if (!face_check || !vert_check)	this->checkFaces<3>();
	if (pvector.size() == 0) this->makePvector<3>(pvector);
	int answer = 0, count = 0;
	for (unsigned int i = 0; i < pvector.size(); ++i) {
		count += pvector[i]; // number of faces
		answer += pvector[i] * i; // number of edges
	}
	return float(answer) / count;
}

template <int dim>
void printCSV(std::ofstream& o, const Grain& g, const int& p, const int& N) {
	if (dim==2) {
		if ( g.isExcluded() ) {
			Point<int> x = g.getCentroid();
			o << g.getID();
			for (int d=0; d<dim; ++d) o << ',' << x[d];
			o << ',' << g.getVolume() << ",0,1,0,0";
			const vector<int>& pvec = g.getPvector<2>();
			assert(pvec.size() > 1);
			for (int i = 2; i < p; ++i) o << ",0";
			const map<set<id_type>, vector<Point<int> > >& edges = g.getEdges();
			for (map<set<id_type>, vector<Point<int> > >::const_iterator itr = edges.begin(); itr != edges.end(); ++itr) o << ',' << *(itr->first.begin());
			for (int i = edges.size(); i < N; ++i) o << ',';
			o << '\n';
		} else if ( g.numVerts() < 1 || g.numEdges() < 2 ) return;
		else {
			Point<int> x = g.getCentroid();
			o << g.getID();
			for (int d=0; d<dim; ++d) o << ',' << x[d];
			o << ",0," << g.getVolume() << ',' << g.numVerts() << ',' << g.numFaces() << ',' << g.numEdges() << ',' << g.getEuler();
			const vector<int>& pvec = g.getPvector<2>();
			assert(pvec.size() > 1);
			for (unsigned int i = 2; i < pvec.size(); ++i) o << ',' << pvec[i];
			for (int i = int(pvec.size()); i < p; ++i) o << ",0";
			const map<set<id_type>, vector<Point<int> > >& edges = g.getEdges();
			for (map<set<id_type>, vector<Point<int> > >::const_iterator itr = edges.begin(); itr != edges.end(); ++itr) o << ',' << *(itr->first.begin());
			for (int i = edges.size(); i < N; ++i) o << ',';
			o << '\n';
		}
	} else {
		if ( g.isExcluded() && g.numFaces() > 0 ) {
			Point<int> x = g.getCentroid();
			o << g.getID();
			for (int d=0; d<dim; ++d) o << ',' << x[d];
			o << ',' << g.getVolume() << ",0," << g.numFaces() << ",0,0";
			const vector<int>& pvec = g.getPvector<3>();
			assert(pvec.size() > 1);
			for (int i = 2; i < p; ++i) o << ",0";
			const set<id_type>& faces = g.getFaces();
			for (set<id_type>::iterator itr = faces.begin(); itr != faces.end(); ++itr) o << ',' << *itr;
			for (int i = faces.size(); i < N; ++i) o << ',';
			o << '\n';
		} else if ( g.numVerts() < 2 || g.numEdges() < 3 || g.numFaces() < 3 ) return; // Fewer F, V, or E than a Brazil nut? Invalid.
		else {
			Point<int> x = g.getCentroid();
			o << g.getID();
			for (int d=0; d<dim; ++d) o << ',' << x[d];
			o << ',' << g.getVolume() << ',' << g.numVerts() << ',' << g.numFaces() << ',' << g.numEdges() << ',' << g.getEuler();
			const vector<int>& pvec = g.getPvector<3>();
			assert(pvec.size() > 1);
			for (unsigned int i = 2; i < pvec.size(); ++i) o << ',' << pvec[i];
			for (int i = int(pvec.size()); i < p; ++i) o << ",0";
			const set<id_type>& faces = g.getFaces();
			for (set<id_type>::iterator itr = faces.begin(); itr != faces.end(); ++itr) o << ',' << *itr;
			for (int i = faces.size(); i < N; ++i) o << ',';
			o << '\n';
		}
	}
}
