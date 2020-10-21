/* graingrowth.cpp
** Algorithms for 2D and 3D isotropic Monte Carlo grain growth
** Ghost communication is performed on 1/4 of all boundaries each time.
** Parallel algorithm: Wright, Steven A., et al. "Potts-model grain growth simulations: Parallel algorithms and applications." SAND Report (1997): 1925.
**
** Questions/comments to tany3@rpi.edu (Yixuan Tan)
*/

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
#include<cmath>
#include<vector>
#include"MMSP.hpp"
#include"graingrowth.hpp"

namespace MMSP
{
void generate(int dim, const char* filename)
{
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.
	if (dim == 1) {
		int L=1024;
		GRID1D initGrid(0, 0, L);

		for (int i = 0; i < nodes(initGrid); i++)
			initGrid(i) = rand() % 100;

		output(initGrid, filename);
	}

	if (dim == 2) {
		int L=256;
		GRID2D initGrid(0, 0, 2*L, 0, L);

		for (int i = 0; i < nodes(initGrid); i++)
			initGrid(i) = rand() % 20;

		output(initGrid, filename);
	} else if (dim == 3) {
		int L=64;
		GRID3D initGrid(0, 0, 2*L, 0, L, 0, L/4);

		for (int i = 0; i < nodes(initGrid); i++)
			initGrid(i) = rand() % 20;

		output(initGrid, filename);
	}
}

template <int dim> bool isOutsideDomain(const grid<dim,int>& mcGrid, const vector<int>& x)
{
	bool outside_domain = false;
	for (int i = 0; i < dim; i++) {
		if (x[i] < x0(mcGrid, i) || x[i] > x1(mcGrid, i)) {
			outside_domain = true;
			break;
		}
	}
	return outside_domain;
}

template <int dim> void update(grid<dim, int>& mcGrid, int steps)
{
	int rank = 0;
	#ifdef MPI_VERSION
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	ghostswap(mcGrid);

	/*---------------generate cells------------------*/
	int dimension_length = 0, num_lattice_cells = 1;
	vector<int> lattice_cells_each_dimension(dim,0);
	for (int i = 0; i < dim; i++) {
		dimension_length = x1(mcGrid, i) - x0(mcGrid, i);
		if (x0(mcGrid, 0) % 2 == 0)
			lattice_cells_each_dimension[i] = dimension_length / 2 + 1;
		else
			lattice_cells_each_dimension[i] = 1 + (dimension_length % 2 == 0 ? dimension_length / 2 : dimension_length / 2 + 1);
		num_lattice_cells *= lattice_cells_each_dimension[i];
	}

	vector<int> x(dim, 0);
	vector<int> x_prim(dim, 0);
	vector<int> initial_coordinates(dim,0);

	vector<int> num_grids_to_flip(int(pow(2, dim)),0);
	vector<int> first_cell_start_coordinates(dim,0);
	for (int kk = 0; kk < dim; kk++)
		first_cell_start_coordinates[kk] = x0(mcGrid, kk);
	for (int i = 0; i < dim; i++)
		if (x0(mcGrid, i) % 2 != 0)
			first_cell_start_coordinates[i]--;


	for (int j = 0; j < num_lattice_cells; j++) {
		vector<int> cell_coords_selected(dim,0);
		if (dim == 2) {
			cell_coords_selected[dim - 1] = j % lattice_cells_each_dimension[dim - 1]; //1-indexed
			cell_coords_selected[0] = (j / lattice_cells_each_dimension[dim - 1]);
		} else if (dim == 3) {
			cell_coords_selected[dim - 1] = j % lattice_cells_each_dimension[dim - 1]; //1-indexed
			cell_coords_selected[1] = (j / lattice_cells_each_dimension[dim - 1]) % lattice_cells_each_dimension[1];
			cell_coords_selected[0] = ( j / lattice_cells_each_dimension[dim - 1] ) / lattice_cells_each_dimension[1];
		}
		for (int i = 0; i < dim; i++)
			x[i] = first_cell_start_coordinates[i] + 2 * cell_coords_selected[i];

		if (dim == 2) {
			x_prim = x;
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[0] += 1;

			x_prim = x;
			x_prim[1] = x[1] + 1; //0,1
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[1] += 1;

			x_prim = x;
			x_prim[0] = x[0] + 1; //1,0
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[2] += 1;

			x_prim = x;
			x_prim[0] = x[0] + 1;
			x_prim[1] = x[1] + 1; //1,1
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[3] += 1;

		} else if (dim == 3) {
			x_prim = x;//0,0,0
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[0] += 1;

			x_prim = x;
			x_prim[2] = x[2] + 1; //0,0,1
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[1] += 1;

			x_prim = x;
			x_prim[1] = x[1] + 1; //0,1,0
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[2] += 1;

			x_prim = x;
			x_prim[2] = x[2] + 1;
			x_prim[1] = x[1] + 1; //0,1,1
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[3] += 1;

			x_prim = x;
			x_prim[0] = x[0] + 1; //1,0,0
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[4] += 1;

			x_prim = x;
			x_prim[2] = x[2] + 1;
			x_prim[0] = x[0] + 1; //1,0,1
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[5] += 1;

			x_prim = x;
			x_prim[1] = x[1] + 1;
			x_prim[0] = x[0] + 1; //1,1,0
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[6] += 1;

			x_prim = x;
			x_prim[2] = x[2] + 1;
			x_prim[1] = x[1] + 1;
			x_prim[0] = x[0] + 1; //1,1,1
			if (!isOutsideDomain<dim>(mcGrid, x_prim)) num_grids_to_flip[7] += 1;
		}
	}// for int j

	for (int k = 0; k < dim; k++)
		initial_coordinates[k] = x0(mcGrid, k);
	for (int i = 0; i < dim; i++)
		if (x0(mcGrid, i) % 2 != 0)
			initial_coordinates[i]--;

	for (int step = 0; step < steps; step++) {
		if (rank == 0)
			print_progress(step, steps);
		int num_sublattices = 0;
		if (dim == 2) num_sublattices = 4;
		else if (dim == 3) num_sublattices = 8;
		for (int sublattice = 0; sublattice < num_sublattices; sublattice++) {

			vector<int> x(dim, 0);
			// This particular algorithm requires that srand() be called here.
			unsigned long seed=time(NULL);
			#ifdef MPI_VERSION
			MPI_Bcast(&seed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
			#endif
			srand(seed); // Also, time(NULL)+rank is an INCORRECT seed for this purpose.

			for (int hh = 0; hh < num_grids_to_flip[sublattice]; hh++) {
				int cell_numbering = rand() % (num_lattice_cells); //choose a cell to flip, from 0 to num_cells_in_thread-1
				vector<int> cell_coords_selected(dim,0);
				if (dim == 2) {
					cell_coords_selected[dim - 1] = cell_numbering % lattice_cells_each_dimension[dim - 1]; //1-indexed
					cell_coords_selected[0] = (cell_numbering / lattice_cells_each_dimension[dim - 1]);
				} else if (dim == 3) {
					cell_coords_selected[dim - 1] = cell_numbering % lattice_cells_each_dimension[dim - 1]; //1-indexed
					cell_coords_selected[1] = (cell_numbering / lattice_cells_each_dimension[dim - 1]) % lattice_cells_each_dimension[1];
					cell_coords_selected[0] = ( cell_numbering / lattice_cells_each_dimension[dim - 1] ) / lattice_cells_each_dimension[1];
				}
				for (int i = 0; i < dim; i++) {
					x[i] = first_cell_start_coordinates[i] + 2 * cell_coords_selected[i];
				}

				if (dim == 2) {
					switch(sublattice) {
					case 0:
						break;// 0,0
					case 1:
						x[1]++;
						break; //0,1
					case 2:
						x[0]++;
						break; //1,0
					case 3:
						x[0]++;
						x[1]++;
						break; //1,1
					}
				} else if (dim == 3) {
					switch(sublattice) {
					case 0:
						break;// 0,0,0
					case 1:
						x[2]++;
						break; //0,0,1
					case 2:
						x[1]++;
						break; //0,1,0
					case 3:
						x[2]++;
						x[1]++;
						break; //0,1,1
					case 4:
						x[0]++;
						break; //1,0,0
					case 5:
						x[2]++;
						x[0]++;
						break; //1,0,1
					case 6:
						x[1]++;
						x[0]++;
						break; //1,1,0
					case 7:
						x[2]++;
						x[1]++;
						x[0]++; //1,1,1
						break;
					}
				}

				bool site_outside_domain = false;
				for (int i = 0; i < dim; i++) {
					if (x[i] < x0(mcGrid, i) || x[i] >= x1(mcGrid, i)) {
						site_outside_domain = true;
						break;//break from the for int i loop
					}
				}
				if (site_outside_domain == true) {
					hh--;
					continue; //continue the int hh loop
				}

				double site_selection_probability = 1.0;

				double rd = double(rand()) / double(RAND_MAX);
				if (rd > site_selection_probability) continue; //this site wont be selected

				int spin1 = mcGrid(x);
				// determine neighboring spins
				vector<int> r(dim, 0);
				std::vector<int> neighbors;
				neighbors.clear();
				unsigned int num_same_neighbors = 0;
				if (dim == 2) {
					for (int i = -1; i <= 1; i++) {
						for (int j = -1; j <= 1; j++) {
							if (!(i == 0 && j == 0)) {
								r[0] = x[0] + i;
								r[1] = x[1] + j;
								if (r[0] < g0(mcGrid, 0) || r[0] >= g1(mcGrid, 0) ||
								    r[1] < g0(mcGrid, 1) || r[1] >= g1(mcGrid, 1) ) // not periodic BC
									continue;// neighbor outside the global boundary, skip it.
								int spin = mcGrid(r);
								neighbors.push_back(spin);
								if (spin == spin1)
									num_same_neighbors++;
							}
						}
					}
				} else if (dim == 3) {
					for (int i = -1; i <= 1; i++) {
						for (int j = -1; j <= 1; j++) {
							for (int k = -1; k <= 1; k++) {
								if (!(i == 0 && j == 0 && k == 0)) {
									r[0] = x[0] + i;
									r[1] = x[1] + j;
									r[2] = x[2] + k;
									if (r[0] < g0(mcGrid, 0) || r[0] >= g1(mcGrid, 0) ||
									    r[1] < g0(mcGrid, 1) || r[1] >= g1(mcGrid, 1) ||
									    r[2] < g0(mcGrid, 2) || r[2] >= g1(mcGrid, 2)) // not periodic BC
										continue;// neighbor outside the global boundary, skip it.
									int spin = mcGrid(r);
									neighbors.push_back(spin);
									if (spin == spin1)
										num_same_neighbors++;
								}
							}
						}
					}
				}

				//check if inside a grain
				if (num_same_neighbors == neighbors.size()) { //inside a grain
					continue;//continue for
				}
				//choose a random neighbor spin
				int spin2 = neighbors[rand() % neighbors.size()];
				// choose a random spin from Q states
				if (spin1 != spin2) {
					// compute energy change
					double dE = 0.0;
					if (dim == 2) {
						for (int i = -1; i <= 1; i++) {
							for (int j = -1; j <= 1; j++) {
								if (!(i == 0 && j == 0)) {
									r[0] = x[0] + i;
									r[1] = x[1] + j;
									if (r[0] < g0(mcGrid, 0) || r[0] >= g1(mcGrid, 0) ||
									    r[1] < g0(mcGrid, 1) || r[1] >= g1(mcGrid, 1) ) // not periodic BC
										continue;// neighbor outside the global boundary, skip it.
									int spin = mcGrid(r);
									dE += 0.5 * ((spin != spin2) - (spin != spin1));
								}// if (!(i==0 && j==0))
							}
						}
					}
					if (dim == 3) {
						for (int i = -1; i <= 1; i++) {
							for (int j = -1; j <= 1; j++) {
								for (int k = -1; k <= 1; k++) {
									if (!(i == 0 && j == 0 && k == 0)) {
										r[0] = x[0] + i;
										r[1] = x[1] + j;
										r[2] = x[2] + k;
										if (r[0] < g0(mcGrid, 0) || r[0] >= g1(mcGrid, 0) ||
										    r[1] < g0(mcGrid, 1) || r[1] >= g1(mcGrid, 1) ||
										    r[2] < g0(mcGrid, 2) || r[2] >= g1(mcGrid, 2)) // not periodic BC
											continue;// neighbor outside the global boundary, skip it.
										int spin = mcGrid(r);
										dE += 0.5 * ((spin != spin2) - (spin != spin1));
									}
								}
							}
						}
					}
					// attempt a spin flip
					double r = double(rand()) / double(RAND_MAX);
					double kT = 1.3803288e-23 * 273.0;

					if (dE <= 0.0)
						mcGrid(x) = spin2;
					else if (r < exp(-dE / kT))
						mcGrid(x) = spin2;
				}//spin1!=spin2
			} // hh

			#ifdef MPI_VERSION
			MPI_Barrier(MPI_COMM_WORLD);
			#endif
			ghostswap(mcGrid, sublattice); // once looped over a "color", ghostswap.
		}//loop over sublattice
	}//loop over step
}//update

}

#endif
#include"MMSP.main.hpp"
