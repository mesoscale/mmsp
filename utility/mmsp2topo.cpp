// File:    mmsp2topo.cpp
// Purpose: reads MMSP data file, analyzes grid for grain topologies
// Output:  Comma-separated-values file summarizing p-vector, vertex connectivity, Euler characteristic

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <utility>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>

#include "point.hpp"
#include "topology.cpp"
#include "MMSP.hpp"

using namespace std;

template <int dim, class T>
Point<T> getPoint(const MMSP::vector<T>& v) {
	Point<T> answer(v[0], v[1], (dim == 3) ? v[2] : 0);
	return answer;
}

// Scan neighboring voxels, adding neighboring grain ID's to a temporary map.
template <int dim, class T>
int check_neighbors(id_type me, const MMSP::vector<int>& v, const MMSP::grid<dim, MMSP::sparse<T> >& data, set<id_type>& adjacent) {
	adjacent.clear();
	MMSP::vector<int> x = v;
	// Check adjacent points for other grain ID's
	for (x[0] = v[0] - 1; x[0] < v[0] + 2; ++x[0]) {
		if (x[0] < MMSP::x0(data, 0) || x[0] >= MMSP::x1(data, 0)) continue;
		for (x[1] = v[1] - 1; x[1] < v[1] + 2; ++x[1]) {
			if (x[1] < MMSP::x0(data, 1) || x[1] >= MMSP::x1(data, 1)) continue;
			if (dim==2) {
				if ( x[0] == v[0] && x[1] == v[1] ) continue;
				else if ( data(x).grain_id() == me ) continue;
				// Add 1st and 2nd nearest neighbors.
				else adjacent.insert( data(x).grain_id() );
			} else {
				for (x[2] = v[2] - 1; x[2] < v[2] + 2; ++x[2]) {
					if (x[2] < MMSP::x0(data, 2) || x[2] >= MMSP::x1(data, 2)) continue;
					if ( x[0] == v[0] && x[1] == v[1] && x[2] == v[2] ) continue;
					else if ( data(x).grain_id() == me ) continue;
					// Skip the 8 corners of this 27-voxel domain. This may skip some valid edge points, but hopefully not all of them.
					else if ( (x[0] == v[0] - 1 || x[0] == v[0] + 1) && (x[1] == v[1] - 1 || x[1] == v[1] + 1) && (x[2] == v[2] - 1 || x[2] == v[2] + 1) ) continue;
					// Add 1st and 2nd nearest neighbors.
					else adjacent.insert( data(x).grain_id() );
				}
			}
		}
	}
	return adjacent.size();
}

int main(int argc, char* argv[]) {
	std::time_t rawtime;
	time_t starttime = std::time( &rawtime );
	time_t timeinfo;
	stringstream timestr;

	if ( argc != 3 ) {
		cout << "Usage: " << argv[0] << " data.dat output.csv\n";
		return ( 1 );
	}

	//  +=========================== Read data file; Define map variables =============================+

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[1] << ".\n\n";
		exit(-1);
	}

	// read data type
	std::string type;
	getline(input, type, '\n');

	// grid type error check
	if (type.substr(0, 4) != "grid") {
		std::cerr << "File input error: " << argv[1] << " does not contain grid data." << std::endl;
		exit(-1);
	} else if (type.substr(5, 6) != "sparse") {
		std::cerr << "File input error: " << argv[1] << " grid does not contain sparse matrices." << std::endl;
		exit(-1);
	} else if (type.substr(12, 5) != "float") {
		std::cerr << "File input error: " << argv[1] << " sparse matrix does not contain single-precision floats." << std::endl;
		exit(-1);
	}

	// read grid dimension
	int dim;
	input >> dim;
	if (dim<2 || dim>3) {
		std::cerr << "ERROR in " << argv[1] << ": "<<dim<<"-dimensional topology is not implemented.\n" << std::endl;
		exit(-1);
	}

	float min_phi = 1., max_phi = 0.;
	map<int, int> topologies; // number of vertices, edges, etc.
	map<int, int> exclusions;
	exclusions[138] = 0;
	exclusions[177] = 0;
	exclusions[202] = 0;
	exclusions[244] = 0;
	exclusions[262] = 0;

	// Read grain ID's from the VTK file. Populate the data matrix.
	map<id_type, Grain> grains;

	// global grid dimensions and spacing
	int* gmin = new int[dim];
	int* gmax = new int[dim];
	float* gscale = new float[dim];

	if (dim==2) {
		// construct grid object
		MMSP::grid<2, MMSP::sparse<float> > grid(argv[1]);
		timeinfo = std::time( &rawtime );
		timestr << "File IO:    " << timeinfo - starttime << " sec. ";
		#ifdef DEBUG
		cout << '\n' << flush;
		#endif
		starttime = timeinfo;

		for (int d=0; d<dim; ++d) {
			gmin[d]=g0(grid,d);
			gmax[d]=g1(grid,d);
			gscale[d]=dx(grid,d);
		}

		// Populate the data matrix from MMSP data.
		for (int i = 0; i < nodes(grid); ++i) {
			MMSP::vector<int> x = position(grid, i);
			id_type major_phase(grid(i).grain_id());
			float phi = grid(i)[major_phase];
			if (phi < min_phi) min_phi = phi;
			else if (phi > max_phi) max_phi = phi;
			map<id_type, Grain>::iterator itr = grains.insert( make_pair(major_phase, Grain(major_phase)) ).first; // second = bool indicating success
			itr->second.addMass(phi, Point<float>(x[0], x[1]));
			int n = length(grid(i));
			for (int h = 0; h < n; ++h) {
				id_type minor_phase(MMSP::index(grid(i), h));
				if (minor_phase == major_phase) continue;
				map<id_type, Grain>::iterator tmp = grains.insert( make_pair(minor_phase, Grain(minor_phase)) ).first;
				float phi = grid(i)[minor_phase];
				tmp->second.addMass(phi, Point<float>(x[0], x[1]));
			}
			set<id_type> neighbors;
			check_neighbors<2,float>(major_phase, x, grid, neighbors);
			++topologies[neighbors.size()];
			if (int(neighbors.size()) == dim-1) {
				// voxel sits on a grain boundary
				itr = grains.find(major_phase);
				itr->second.updateEdges(neighbors, Point<int>(x[0], x[1]));
			}	else if (int(neighbors.size()) == dim) {
				// voxel sits on a grain boundary vertex
				for (set<id_type>::const_iterator j=neighbors.begin(); j!=neighbors.end(); ++j) {
					itr = grains.find(*j);
					itr->second.updateVertices(neighbors);
				}
			}
		}
		#ifdef DEBUG
			cout << "\n Neighbors  Count\n";
			for (map<int, int>::iterator itr = topologies.begin(); itr != topologies.end(); ++itr) cout << " " << setw(9) << right << itr->first << "  " << itr->second << '\n';
			cout << '\n';
		#endif
	} else if (dim==3) {
		// construct grid object
		MMSP::grid<3, MMSP::sparse<float> > grid(argv[1]);
		timeinfo = std::time( &rawtime );
		timestr << "File IO:    " << timeinfo - starttime << " sec. ";
		#ifdef DEBUG
		cout << '\n' << flush;
		#endif
		starttime = timeinfo;

		for (int d=0; d<dim; ++d) {
			gmin[d]=g0(grid,d);
			gmax[d]=g1(grid,d);
			gscale[d]=dx(grid,d);
		}

		// Populate the data matrix from MMSP data.
		for (int i = 0; i < nodes(grid); ++i) {
			MMSP::vector<int> x = position(grid, i);
			id_type major_phase(grid(i).grain_id());
			float phi = grid(i)[major_phase];
			if (phi < min_phi) min_phi = phi;
			else if (phi > max_phi) max_phi = phi;
			map<id_type, Grain>::iterator itr = grains.insert( make_pair(major_phase, Grain(major_phase)) ).first; // second = bool indicating success
			itr->second.addMass(phi, Point<float>(x[0], x[1], x[2]));
			int n = length(grid(i));
			for (int h = 0; h < n; ++h) {
				id_type minor_phase(MMSP::index(grid(i), h));
				if (minor_phase == major_phase) continue;
				map<id_type, Grain>::iterator tmp = grains.insert( make_pair(minor_phase, Grain(minor_phase)) ).first;
				float phi = grid(i)[minor_phase];
				tmp->second.addMass(phi, Point<float>(x[0], x[1], x[2]));
			}
			set<id_type> neighbors;
			check_neighbors<3,float>(major_phase, x, grid, neighbors);
			++topologies[neighbors.size()];
			if (int(neighbors.size()) == dim-1) {
				// voxel sits on a triple-line between grains
				itr = grains.find(major_phase);
				itr->second.updateEdges(neighbors, Point<int>(x[0], x[1], x[2]));
				if (dim==3)
					for (set<id_type>::iterator sitr = neighbors.begin(); sitr != neighbors.end(); ++sitr)
						itr->second.addNeighbor(*sitr);
			}
		}
		#ifdef DEBUG
			cout << "\n Neighbors  Count\n";
			for (map<int, int>::iterator itr = topologies.begin(); itr != topologies.end(); ++itr) cout << " " << setw(9) << right << itr->first << "  " << itr->second << '\n';
			cout << '\n';
		#endif
	}


	//  +=========================== Distinguish Bulk from Surface grains ===============================+

	float average_volume = 0.;
	int BulkGrains = 0;
	for (map<id_type, Grain>::iterator itr = grains.begin(); itr != grains.end(); ++itr) {
		if ( itr->second.isExcluded() ) continue;
		float vol = itr->second.getVolume();
		if (vol < 1.0) {
			if (!(itr->second.isExcluded())) {
				++exclusions[177];
				itr->second.setExclusion();
			}
			continue;
		}
		++BulkGrains;
		if (dim==2) average_volume += gscale[0] * gscale[1] * vol; // square microns
		else average_volume += gscale[0] * gscale[1] * gscale[2] * vol; // cubic microns
	}
	average_volume /= BulkGrains;
	float average_radius;
	if (dim==2) {
		average_radius= sqrt(average_volume/M_PI); // Circle-equivalent grain radius
	} else {
		average_radius= pow(3 * average_volume / (4 * M_PI), 1. / 3.); // Sphere-equivalent grain radius
	}
	Point<float> average_radii(average_radius / gscale[0], average_radius / gscale[1], 0.0f);
	if (dim==3) average_radii.z=average_radius / gscale[2]; // voxel units

	Point<float> centroid(0, 0, 0); // should be (L/2, L/2, L/2)... right? ... and so it is, roughly. Depends on exclusion.
	Point<int> min(1000, 1000);
	if (dim==3) min.z=1000;
	Point<int> max(0, 0, 0);
	BulkGrains = 0;
	for (map<id_type, Grain>::iterator itr = grains.begin(); itr != grains.end(); ++itr) {
		Point<float> center = itr->second.getCentroid();
		if (itr->second.isExcluded()) continue;
		++BulkGrains;
		centroid += center;
		for (int d=0; d<dim; ++d) {
			if (center[d] < min[d]) min[d] = center[d];
			else if (center[d] > max[d]) max[d] = center[d];
		}
	}
	centroid /= BulkGrains;
	cout << argv[1] << ": " << BulkGrains << " bulk grains, " << grains.size() - BulkGrains << " surface grains. Phase spanned [" << min_phi << ", " << max_phi << "].\n";
	cout << "The average grain has volume " << setprecision(0) << fixed << average_volume << " cubic microns, sphere-equivalent radius " << average_radius << " units, centered at " << setprecision(0) << fixed << centroid << ".\n";
	cout << "Centers-of-Mass span " << setprecision(0) << fixed << min << " -- " << setprecision(0) << fixed << max << ".\n";
	timeinfo = std::time( &rawtime );
	timestr << "Setup:    " << timeinfo - starttime << " sec. ";
	#ifdef DEBUG
	cout << '\n' << flush;
	#endif
	starttime = timeinfo;

	delete [] gmin; gmin=NULL;
	delete [] gmax; gmax=NULL;
	delete [] gscale; gscale=NULL;

	//  +=========================== Process the data matrix, mapping each Edge voxel to a grain =============================+

	int maxp = 0, maxN = 0;

	#ifdef DEBUG
	cout << "\nAnalyzing grains.\n";
	#endif

	for(map<id_type, Grain>::iterator gitr = grains.begin(); gitr != grains.end(); ++gitr) {
		// Sanitize topologies
		#ifdef DEBUG
		cout << '\n' << setw(5) << right << gitr->first << ": " << "a" << flush;
		#endif

		// Infer vertices from edges
		std::vector<int> stats;
		if (dim==2) stats = gitr->second.computeTopology<2>();
		else stats = gitr->second.computeTopology<3>();
		assert(stats.size() == 4);
		if (stats[0] == -1) {
			if (!(gitr->second.isExcluded())) {
				++exclusions[244];
				gitr->second.setExclusion();
			}
			#ifdef DEBUG
			cout << " No faces: excluding. " << flush;
			#endif
		} else if (stats[0] != 0) {
			#ifdef DEBUG
			cout << ' ' << stats[0] << " improper faces were erased, " << gitr->second.numFaces() << " remain." << flush;
			#endif
		}
		if (stats[2] > maxp) maxp = stats[2]; // pvector.size() == maxp
		if (stats[3] > maxN) maxN = stats[3]; // neighbors.size() == maxN
		#ifdef DEBUG
		cout << "b" << flush;
		#endif
		if (gitr->second.numEdges() == 0) {
			if (!(gitr->second.isExcluded())) {
				++exclusions[262];
				gitr->second.setExclusion();
			}
#ifdef DEBUG
			cout << " No edges: Excluding." << flush;
#endif
			continue;
		}
#ifdef DEBUG
		cout << "c" << flush;
		cout << " x=" << gitr->second.getEuler();
#endif
	}
	timeinfo = std::time( &rawtime );
	timestr << "Analysis: " << timeinfo - starttime << " sec.\n";

#ifdef DEBUG
	cout << "\n\nWriting output." << endl;
#endif
	// At this point, the global maximum p-vector and neighbor lists are known.
	unsigned int euler_error = 0, smith_error = 0, plateau_error = 0;
	ofstream ofile( argv[2] );
	ofile << "Grain ID,Centroid(X),Centroid(Y),Centroid(Z),Bulk Volume,Vertices,Faces,Edges,Euler";
	for (int i = 2; i < maxp; ++i) ofile << ",p" << i;
	for (int i = 0; i < maxN; ++i) ofile << ",N" << i;
	ofile << '\n';
	for (map<id_type, Grain>::iterator itr = grains.begin(); itr != grains.end(); ++itr) {
		if (itr->second.isExcluded()) continue;
		if (dim==2) {
			if (int(1) != itr->second.getEuler()) ++euler_error;
			if (!(itr->second.isPlateaunic<2>())) ++plateau_error;
			printCSV<2>(ofile, itr->second, maxp, maxN); // finishes with a newline
		} else {
			if (int(2) != (itr->second.getEuler())) ++euler_error;
			if (fabs((6 - 12. / itr->second.numFaces()) - itr->second.getPbar<3>()) > 0.01) ++smith_error;
			if (!(itr->second.isPlateaunic<3>())) ++plateau_error;
			printCSV<3>(ofile, itr->second, maxp, maxN); // finishes with a newline
		}
	}
	ofile << '\n';
	// Print excluded grains
	for (map<id_type, Grain>::iterator itr = grains.begin(); itr != grains.end(); ++itr) {
		if (!(itr->second.isExcluded())) continue;
		if (dim==2) {
			itr->second.estimateCentroid<2>(); // unify bifurcated set
			printCSV<2>(ofile, itr->second, maxp, maxN); // finishes with a newline
		} else {
			itr->second.estimateCentroid<3>(); // unify bifurcated set
			printCSV<3>(ofile, itr->second, maxp, maxN); // finishes with a newline
		}
	}
	ofile << '\n';

	//cout<<setw(4)<<right<<euler_error<<" grains of "<<setw(4)<<right<<grains.size()<<" violated the Euler characteristic (inference from edges).\n";
	cout << setw(4) << right << smith_error << " grains of " << setw(4) << right << BulkGrains << " differed from Smith's formula for <p>.\n";
	cout << setw(4) << right << plateau_error << " grains of " << setw(4) << right << grains.size() << " violated Plateau's Law (3-regularity).\n";
	//cout<<"Exclusions due to 138: "<<exclusions[138]<<", 177: "<<exclusions[177]<<", 202: "<<exclusions[202]<<", 244: "<<exclusions[244]<<", 262: "<<exclusions[262]<<".\n";
	cout << "Exclusions due to Domain-boundary: " << exclusions[138]
	     << ", Zero-volume: " << exclusions[177]
	     << ", Exclusion-zone: " << exclusions[202]
	     << ", Zero-faces: " << exclusions[244]
	     << ", Zero-edges: " << exclusions[262]
	     << '.' << endl;

	cout << timestr.rdbuf() << '\n';
	return 0;
}
