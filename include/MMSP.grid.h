// MMSP.grid.h
// MMSP grid class definition
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MMSP_GRID
#define MMSP_GRID

#include "MMSP.utility.h"
#include "MMSP.scalar.h"
#include "MMSP.vector.h"
#include "MMSP.sparse.h"


namespace MMSP
{


// declaration of grid class
template <int dim, typename T> class grid;

// grid utility functions
template <int dim, typename T> int nodes(const grid<dim, T>& GRID)
{
	return nodes(GRID);
}
template <int dim, typename T> int fields(const grid<dim, T>& GRID)
{
	return fields(GRID);
}
template <int dim, typename T> int ghosts(const grid<dim, T>& GRID)
{
	return ghosts(GRID);
}
template <int dim, typename T> int g0(const grid<dim, T>& GRID, int i)
{
	return g0(GRID, i);
}
template <int dim, typename T> int g1(const grid<dim, T>& GRID, int i)
{
	return g1(GRID, i);
}
template <int dim, typename T> int b0(const grid<dim, T>& GRID, int i)
{
	return b0(GRID, i);
}
template <int dim, typename T> int b1(const grid<dim, T>& GRID, int i)
{
	return b1(GRID, i);
}
template <int dim, typename T> int& b0(grid<dim, T>& GRID, int i)
{
	return b0(GRID, i);
}
template <int dim, typename T> int& b1(grid<dim, T>& GRID, int i)
{
	return b1(GRID, i);
}

// grid utility functions (all directions)
template <int dim, typename T> int x0(const grid<dim, T>& GRID, int i)
{
	return x0(GRID, i);
}
template <int dim, typename T> int x1(const grid<dim, T>& GRID, int i)
{
	return x1(GRID, i);
}
template <int dim, typename T> int xmin(const grid<dim, T>& GRID, int i)
{
	return xmin(GRID, i);
}
template <int dim, typename T> int xmax(const grid<dim, T>& GRID, int i)
{
	return xmax(GRID, i);
}
template <int dim, typename T> int xlength(const grid<dim, T>& GRID, int i)
{
	return xlength(GRID, i);
}
template <int dim, typename T> double dx(const grid<dim, T>& GRID, int i)
{
	return dx(GRID, i);
}
template <int dim, typename T> double& dx(grid<dim, T>& GRID, int i)
{
	return dx(GRID, i);
}

// grid utility functions (x direction)
template <int dim, typename T> int x0(const grid<dim, T>& GRID)
{
	return x0(GRID);
}
template <int dim, typename T> int x1(const grid<dim, T>& GRID)
{
	return x1(GRID);
}
template <int dim, typename T> int xmin(const grid<dim, T>& GRID)
{
	return xmin(GRID);
}
template <int dim, typename T> int xmax(const grid<dim, T>& GRID)
{
	return xmax(GRID);
}
template <int dim, typename T> int xlength(const grid<dim, T>& GRID)
{
	return xlength(GRID);
}
template <int dim, typename T> double dx(const grid<dim, T>& GRID)
{
	return dx(GRID);
}
template <int dim, typename T> double& dx(grid<dim, T>& GRID)
{
	return dx(GRID);
}

// grid utility functions (y direction)
template <int dim, typename T> int y0(const grid<dim, T>& GRID)
{
	return y0(GRID);
}
template <int dim, typename T> int y1(const grid<dim, T>& GRID)
{
	return y1(GRID);
}
template <int dim, typename T> int ymin(const grid<dim, T>& GRID)
{
	return ymin(GRID);
}
template <int dim, typename T> int ymax(const grid<dim, T>& GRID)
{
	return ymax(GRID);
}
template <int dim, typename T> int ylength(const grid<dim, T>& GRID)
{
	return ylength(GRID);
}
template <int dim, typename T> double dy(const grid<dim, T>& GRID)
{
	return dy(GRID);
}
template <int dim, typename T> double& dy(grid<dim, T>& GRID)
{
	return dy(GRID);
}

// grid utility functions (z direction)
template <int dim, typename T> int z0(const grid<dim, T>& GRID)
{
	return z0(GRID);
}
template <int dim, typename T> int z1(const grid<dim, T>& GRID)
{
	return z1(GRID);
}
template <int dim, typename T> int zmin(const grid<dim, T>& GRID)
{
	return zmin(GRID);
}
template <int dim, typename T> int zmax(const grid<dim, T>& GRID)
{
	return zmax(GRID);
}
template <int dim, typename T> int zlength(const grid<dim, T>& GRID)
{
	return zlength(GRID);
}
template <int dim, typename T> double dz(const grid<dim, T>& GRID)
{
	return dz(GRID);
}
template <int dim, typename T> double& dz(grid<dim, T>& GRID)
{
	return dz(GRID);
}


// instantiation of grid class
template <int dim, typename T>
class grid
{
public:
	// constructors
	grid(int FIELDS, int min[dim], int max[dim], int GHOSTS=1, bool SINGLE=false);

	grid(int FIELDS, ...);

	grid(const grid& GRID);

	template <typename U>
	grid(const grid<dim, U>& GRID, int FIELDS);

	grid(const char* filename, int GHOSTS = 1);

	void setup(bool SINGLE = false);
	// destructor
	~grid();

	// assignment operators

	template<typename U> grid& operator=(const U& value);

	template <typename U> grid& operator=(const grid<dim, U>& GRID);

	template <typename U> grid& operator+=(const grid<dim, U>& GRID);

	template <typename U> grid& operator-=(const grid<dim, U>& GRID);

	// subscript operators
	target < dim - 1, 0, T > operator [](int x) const;

	T& operator()(MMSP::vector<int> x) const;

	T& operator()(int i) const;


	// position utility function
	MMSP::vector<int> position(int i) const;


	// parallelization
	void ghostswap();

	void ghostswap(const int sublattice);

	unsigned long buffer_size_save(const int min[dim], const int max[dim]) const;

	unsigned long buffer_size_save(T* p, int i, const int min[dim], const int max[dim]) const;

	unsigned long to_buffer_save(char* buffer, const int min[dim], const int max[dim]) const;

	unsigned long to_buffer_save(char* buffer, T* p, int i, const int min[dim], const int max[dim]) const;

	unsigned long from_buffer_save(char* buffer, const int min[dim], const int max[dim]);

	unsigned long from_buffer_save(char* buffer, T* p, int i, const int min[dim], const int max[dim]);

	unsigned long buffer_size() const;

	unsigned long buffer_size(const int min[dim], const int max[dim]) const;

	unsigned long buffer_size(T* p, int i, const int min[dim], const int max[dim]) const;

	unsigned long to_buffer(char* buffer) const;

	unsigned long to_buffer(char* buffer, const int min[dim], const int max[dim]) const;

	unsigned long to_buffer(char* buffer, T* p, int i, const int min[dim], const int max[dim]) const;

	unsigned long from_buffer(char* buffer);

	unsigned long from_buffer(char* buffer, const int min[dim], const int max[dim]);

	unsigned long from_buffer(char* buffer, T* p, int i, const int min[dim], const int max[dim]);

	// file I/O
	void input(const char* filename, int GHOSTS = 1, int SINGLE = false);

	void read(std::ifstream& file, int GHOSTS = 1);

	void output(const char* filename) const;

	unsigned long write_buffer(char*& buf) const;

	// grid utility functions
	friend int nodes(const grid& GRID)
	{
		return GRID.nodes;
	}
	friend int fields(const grid& GRID)
	{
		return GRID.fields;
	}
	friend int ghosts(const grid& GRID)
	{
		return GRID.ghosts;
	}
	friend int g0(const grid& GRID, int i)
	{
		return GRID.g0[i];
	}
	friend int g1(const grid& GRID, int i)
	{
		return GRID.g1[i];
	}
	friend int b0(const grid& GRID, int i)
	{
		return GRID.b0[i];
	}
	friend int b1(const grid& GRID, int i)
	{
		return GRID.b1[i];
	}
	friend int& b0(grid& GRID, int i)
	{
		return GRID.b0[i];
	}
	friend int& b1(grid& GRID, int i)
	{
		return GRID.b1[i];
	}

	// grid utility functions (all directions)
	friend int x0(const grid& GRID, int i)
	{
		return GRID.x0[i];
	}
	friend int x1(const grid& GRID, int i)
	{
		return GRID.x1[i];
	}
	friend int xmin(const grid& GRID, int i)
	{
		return GRID.x0[i];
	}
	friend int xmax(const grid& GRID, int i)
	{
		return GRID.x1[i];
	}
	friend int xlength(const grid& GRID, int i)
	{
		return GRID.x1[i] - GRID.x0[i];
	}
	friend double dx(const grid& GRID, int i)
	{
		return GRID.dx[i];
	}
	friend double& dx(grid& GRID, int i)
	{
		return GRID.dx[i];
	}
	friend int N0(const grid& GRID, int i)
	{
		return GRID.n0[i];
	}
	friend int N1(const grid& GRID, int i)
	{
		return GRID.n1[i];
	}
	friend int P0(const grid& GRID, int i)
	{
		return GRID.p0[i];
	}
	friend int P1(const grid& GRID, int i)
	{
		return GRID.p1[i];
	}
	friend int sp(const grid& GRID, int i)
	{
		return GRID.sp[i];
	}

	// grid utility functions (x direction)
	friend int x0(const grid& GRID)
	{
		return GRID.x0[0];
	}
	friend int x1(const grid& GRID)
	{
		return GRID.x1[0];
	}
	friend int xmin(const grid& GRID)
	{
		return GRID.x0[0];
	}
	friend int xmax(const grid& GRID)
	{
		return GRID.x1[0];
	}
	friend int xlength(const grid& GRID)
	{
		return GRID.x1[0] - GRID.x0[0];
	}
	friend double dx(const grid& GRID)
	{
		return GRID.dx[0];
	}
	friend double& dx(grid& GRID)
	{
		return GRID.dx[0];
	}

	// grid utility functions (y direction)
	friend int y0(const grid& GRID)
	{
		return GRID.x0[1];
	}
	friend int y1(const grid& GRID)
	{
		return GRID.x1[1];
	}
	friend int ymin(const grid& GRID)
	{
		return GRID.x0[1];
	}
	friend int ymax(const grid& GRID)
	{
		return GRID.x1[1];
	}
	friend int ylength(const grid& GRID)
	{
		return GRID.x1[1] - GRID.x0[1];
	}
	friend double dy(const grid& GRID)
	{
		return GRID.dx[1];
	}
	friend double& dy(grid& GRID)
	{
		return GRID.dx[1];
	}

	// grid utility functions (z direction)
	friend int z0(const grid& GRID)
	{
		return GRID.x0[2];
	}
	friend int z1(const grid& GRID)
	{
		return GRID.x1[2];
	}
	friend int zmin(const grid& GRID)
	{
		return GRID.x0[2];
	}
	friend int zmax(const grid& GRID)
	{
		return GRID.x1[2];
	}
	friend int zlength(const grid& GRID)
	{
		return GRID.x1[2] - GRID.x0[2];
	}
	friend double dz(const grid& GRID)
	{
		return GRID.dx[2];
	}
	friend double& dz(grid& GRID)
	{
		return GRID.dx[2];
	}


	// utility functions
	void swap(grid& GRID);

	void copy(const grid& GRID);

protected:
	T* data;        // local grid data

	int nodes;			// number of nodes (excluding ghosts)
	int cells;			// number of nodes (including ghosts)
	int fields;			// number of fields
	int ghosts;			// ghost cell depth

	int g0[dim];    // global lower coordinate limit (excluding ghosts)
	int g1[dim];    // global upper coordinate limit (excluding ghosts)

	int x0[dim];    // local lower coordinate limit (excluding ghosts)
	int x1[dim];    // local upper coordinate limit (excluding ghosts)
	int xx[dim];    // local cells in slice (excluding ghosts)

	int s0[dim];    // local lower coordinate limit (including ghosts)
	int s1[dim];    // local upper coordinate limit (including ghosts)
	int sx[dim];    // local cells in slice (including ghosts)

	int b0[dim];    // boundary condition at x0
	int b1[dim];    // boundary condition at x1

	double dx[dim]; // global cell spacing

	int p0[dim];
	int p1[dim];
	int sp[dim];    // global processors in slice

	int n0[dim];    // neighbor processor at x0
	int n1[dim];    // neighbor processor at x1
};


// math functions
template <int dim, typename T> T laplacian(const grid<dim, T>& GRID, const vector<int>& x);

template <int dim, typename T> vector<T> laplacian(const grid<dim, vector<T> >& GRID, const vector<int>& x);

template<int dim, typename T> T laplacian(const grid<dim,vector<T> >& GRID, const vector<int>& x, const int field);

template <int dim, typename T> sparse<T> laplacian(const grid<dim, sparse<T> >& GRID, const vector<int>& x);

template <int dim, typename T> T laplacian(const grid<dim, T>& GRID, int i);

template <int dim, typename T> vector<T> laplacian(const grid<dim, vector<T> >& GRID, int i);


template <int dim, typename T> T laplacian(const grid<dim, vector<T> >& GRID, int i, int f);

template <int dim, typename T> sparse<T> laplacian(const grid<dim, sparse<T> >& GRID, int i);

template <int dim, typename T> vector<T> gradient(const grid<dim, T>& GRID, const vector<int>& x);

template <int dim, typename T> vector<T> gradient(const grid<dim,vector<T> >& GRID, const vector<int>& x, const int field);

template <int dim, typename T> vector<T> grad(const grid<dim, T>& GRID, const vector<int>& x);

template <int dim, typename T> T divergence(const grid<dim, T>& GRID, const vector<int>& x);

template <int dim, typename T> vector<T> divergence(const grid<dim, vector<T> >& GRID, const vector<int>& x);

template <int dim, typename T> T div(const grid<dim, T>& GRID, const vector<int>& x)
{
	return divergence(GRID, x);
}

template <int dim, typename T> vector<T> div(const grid<dim, vector<T> >& GRID, const vector<int>& x)
{
	return divergence(GRID, x);
}

// position utility function
template <int dim, typename T> MMSP::vector<int> position(const grid<dim, T>& GRID, int index)
{
	return GRID.position(index);
}

// parallelization
template <int dim, typename T> void ghostswap(grid<dim, T>& GRID)
{
	GRID.ghostswap();
}

template <int dim, typename T> void ghostswap(grid<dim, T>& GRID, const int sublattice)
{
	GRID.ghostswap(sublattice);
}
// buffer I/O functions
template <int dim, typename T> unsigned long buffer_size(const grid<dim, T>& GRID)
{
	return GRID.buffer_size();
}
template <int dim, typename T> unsigned long to_buffer(const grid<dim, T>& GRID, char* buffer)
{
	return GRID.to_buffer(buffer);
}
template <int dim, typename T> unsigned long from_buffer(grid<dim, T>& GRID, const char* buffer)
{
	return GRID.from_buffer(buffer);
}

// file I/O functions
template <int dim, typename T> void read(grid<dim, T>& GRID, std::ifstream& file)
{
	GRID.read(file);
}
template <int dim, typename T> void write(const grid<dim, T>& GRID, std::ifstream& file)
{
	GRID.write(file);
}
template <int dim, typename T> void input(grid<dim, T>& GRID, const char* filename, int GHOSTS = 1, int SINGLE = false)
{
	GRID.input(filename, GHOSTS, SINGLE);
}
template <int dim, typename T> void output(const grid<dim, T>& GRID, const char* filename)
{
	GRID.output(filename);
}
template <int dim, typename T> unsigned long write_buffer(const grid<dim, T>& GRID, char*& buf)
{
	return GRID.write_buffer(buf);
}

// utility functions
template <int dim, typename T> int length(const grid<dim, T>& GRID)
{
	return nodes(GRID);
}
template <int dim, typename T> void resize(int n, grid<dim, T>& GRID) {}
template <int dim, typename T> void swap(grid<dim, T>& GRID1, grid<dim, T>& GRID2)
{
	GRID1.swap(GRID2);
}
template <int dim, typename T> void copy(grid<dim, T>& GRID1, grid<dim, T>& GRID2)
{
	GRID2.copy(GRID1);
}
template <int dim, typename T> std::string name(const grid<dim, T>& GRID)
{
	return std::string("grid:") + name(T());
}


} // namespace MMSP


#include "MMSP.grid.cpp"

#endif
