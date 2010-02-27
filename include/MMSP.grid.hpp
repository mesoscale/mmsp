// MMSP.grid.hpp
// MMSP grid class definition and implementation 
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MMSP_GRID
#define MMSP_GRID
#include<iostream>
#include<cstdlib>
#include<cstdarg>
#include<sstream>
#include<cmath>

#include"MMSP.utility.hpp"
#include"MMSP.scalar.hpp"
#include"MMSP.vector.hpp"
#include"MMSP.sparse.hpp"

namespace MMSP{

// grid class
template <int dim, typename T>
class grid{ 
public:
	// constructors
	grid(int FIELDS, int min[dim], int max[dim], int GHOSTS=1, int SINGLE=false)
	{
		// set number of fields
		fields = FIELDS;

		// read function arguments
		for (int i=0; i<dim; i++) {
			g0[i] = min[i];
			g1[i] = max[i];
		}

		// set number of ghosts
		ghosts = 0;

		#ifdef MPI_VERSION
		ghosts = GHOSTS;
		#endif

		// setup grid properties
		setup(SINGLE);
	}

	grid(int FIELDS, ...)
	{
		// set number of fields
		fields = FIELDS;

		// read function arguments
		va_list list;
		va_start(list,FIELDS);
		for (int i=0; i<dim; i++) {
			g0[i] = va_arg(list,int);
			g1[i] = va_arg(list,int);
		}

		va_end(list);

		// set number of ghosts
		ghosts = 0;

		#ifdef MPI_VERSION
		ghosts = 1;
		#endif

		// setup grid properties
		setup();
	}

	grid(const grid& GRID)
	{
		// set number of fields
		fields = GRID.fields;

		// read function arguments
		for (int i=0; i<dim; i++) {
			g0[i] = GRID.g0[i];
			g1[i] = GRID.g1[i];
		}

		// set number of ghosts
		ghosts = 0;

		#ifdef MPI_VERSION
		ghosts = GRID.ghosts;
		#endif

		// setup grid properties
		setup();

		// replace defaults
		for (int i=0; i<dim; i++) {
			b0[i] = GRID.b0[i];
			b1[i] = GRID.b1[i];
			dx[i] = GRID.dx[i];
		}
	}

	grid(const grid& GRID, int FIELDS)
	{
		// set number of fields
		fields = FIELDS;

		// read function arguments
		for (int i=0; i<dim; i++) {
			g0[i] = GRID.g0[i];
			g1[i] = GRID.g1[i];
		}

		// set number of ghosts
		ghosts = 0;

		#ifdef MPI_VERSION
		ghosts = GRID.ghosts;
		#endif

		// setup grid properties
		setup();

		// replace defaults
		for (int i=0; i<dim; i++) {
			b0[i] = GRID.b0[i];
			b1[i] = GRID.b1[i];
			dx[i] = GRID.dx[i];
		}
	}

	grid(const char* filename, int GHOSTS=1, int SINGLE=false)
	{
		// initialize data
		data=NULL;

		// read data from file
		input(filename,GHOSTS,SINGLE);

		// ghostswap if necessary
		if (not SINGLE) ghostswap();
	}

	void setup(bool SINGLE=false)
	{
		// setup default grid parameters
		for (int i=0; i<dim; i++) {
			x0[i] = g0[i];
			x1[i] = g1[i];

			b0[i] = periodic;
			b1[i] = periodic;

			dx[i] = 1.0;

			p0[i] = 0;
			p1[i] = 1;
			sp[i] = 1;

			n0[i] = 0;
			n1[i] = 0;
		}

		#ifdef MPI_VERSION
		// get global grid data, set neighbor processors
		int id = MPI::COMM_WORLD.Get_rank();
		int np = MPI::COMM_WORLD.Get_size();

		// if bool SINGLE is set to true,
		// we generate grid data only on proc 0
		if (not SINGLE) {
			// compute integral factors of "np"
			int nfac = 0;
			for (int i=1; i<=np; i++)
				if ((np/i)*i==np) nfac += 1;

			int* factors = new int[nfac];
			nfac = 0;
			for (int i=1; i<=np; i++)
				if ((np/i)*i==np) {
					factors[nfac] = i;
					nfac += 1;
				}

			// compute global slice areas 
			int area[dim];
			for (int i=0; i<dim; i++) {
				area[i] = 1;
				for (int j=0; j<dim; j++)
					if (i!=j) area[i] *= (g1[j]-g0[j]);
			}

			// initialize optimal ghost area
			int minarea = -1;

			// compute all combinations of "dim" factors
			int ncom = 1;
			for (int i=0; i<dim; i++)
				ncom *= nfac;

			for (int i=0; i<ncom; i++) {
				int combo[dim];
				int total = i;
				for (int j=0; j<dim; j++) {
					int slice = 1;
					for (int k=j+1; k<dim; k++)
						slice *= nfac;
					combo[j] = total/slice;
					total -= combo[j]*slice;
				}

				// compute the product of "dim" factors
				int product = 1;
				for (int j=0; j<dim; j++)
					product *= factors[combo[j]];

				// if product equals "np", compute ghost area
				if (product==np) {
					int test = 0;
					for (int j=0; j<dim; j++)
						test += area[j]*(1+factors[combo[j]]);
					// choose optimal (smallest) ghost area
					if (test<minarea or minarea<0) {
						minarea = test;
						for (int k=0; k<dim; k++)
							p1[k] = factors[combo[k]];
					}
				}
			}

			// clean up
			delete [] factors;

			// compute slice sizes
			for (int i=0; i<dim; i++) {
				sp[i] = 1;
				for (int j=i+1; j<dim; j++)
					sp[i] *= p1[j];
			}

			// determine local grid limits
			int pos[dim];

			int total = id;
			for (int i=0; i<dim; i++) {
				// compute position 
				pos[i] = total/sp[i];
				total -= pos[i]*sp[i];

				int nom = (g1[i]-g0[i])/p1[i];
				int rem = (g1[i]-g0[i])-p1[i]*nom;

				// set limits to "nominal" values
				x0[i] = g0[i]+pos[i]*nom;
				x1[i] = x0[i]+nom;

				// adjust limits according to "remainder"
				if (rem>0) {
					if (pos[i]<rem) {
						x0[i] += pos[i];
						x1[i] = x0[i]+nom+1;
					}
					else if (pos[i]==rem) {
						x0[i] += pos[i];
						x1[i] = x0[i]+nom;
					}
					else {
						x0[i] += rem;
						x1[i] = x0[i]+nom;
					}
				}
			}

			// determine neighbor processors
			for (int i=0; i<dim; i++) {
				int npos[dim];
				for (int j=0; j<dim; j++)
					npos[j] = pos[j];

				// set neighbor below
				n0[i] = 0;
				npos[i] = (pos[i]-1+p1[i])%p1[i];
				for (int j=0; j<dim; j++)
					n0[i] += sp[j]*npos[j];

				// set neighbor above
				n1[i] = 0;
				npos[i] = (pos[i]+1)%p1[i];
				for (int j=0; j<dim; j++)
					n1[i] += sp[j]*npos[j];
			}

			// adjust boundary conditions
			for (int i=0; i<dim; i++) {
				if (x0[i]!=g0[i]) b0[i] = parallel;
				if (x1[i]!=g1[i]) b1[i] = parallel;
			}
		}
		#endif

		// compute slice sizes
		for (int i=0; i<dim; i++) {
			s0[i] = x0[i]-ghosts;
			s1[i] = x1[i]+ghosts;
		}
		for (int i=0; i<dim; i++) {
			sx[i] = 1;
			for (int j=i+1; j<dim; j++)
				sx[i] *= (s1[j]-s0[j]);
		}

		// compute number of cells
		cells = sx[0]*(s1[0]-s0[0]);

		// resize data structures 
		data = new T[cells];
		for (int i=0; i<cells; i++)
			resize(data[i],fields);

		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	}


	// destructor
	~grid() {delete [] data;}


	// assignment operators
	grid& operator=(const grid& GRID)
	{
		for (int i=0; i<cells; i++)
			data[i] = GRID.data[i];
	}

	grid& operator=(const T& value)
	{
		for (int i=0; i<cells; i++)
			data[i] = value;
	}

	grid& operator+=(const grid& GRID2)
	{
		for (int i=0; i<cells; i++)
			data[i] += GRID2.data[i];
	}

	grid& operator-=(const grid& GRID2)
	{
		for (int i=0; i<cells; i++)
			data[i] -= GRID2.data[i];
	}


	// subscript operators
	target<dim-1,0,T> operator [](int x)
	{
		check_boundary(x,x0[0],x1[0],b0[0],b1[0]);
		return target<dim-1,0,T>(data+(x-s0[0])*sx[0],s0,sx,x0,x1,b0,b1);
	}

	target<dim-1,0,T> operator [](int x) const
	{
		check_boundary(x,x0[0],x1[0],b0[0],b1[0]);
		return target<dim-1,0,T>(data+(x-s0[0])*sx[0],s0,sx,x0,x1,b0,b1);
	}

	T& operator()(int x[dim]) const
	{
		T* p = data;
		for (int i=0; i<dim; i++) {
			check_boundary(x[i],x0[i],x1[i],b0[i],b1[i]);
			p += (x[i]-s0[i])*sx[i];
		}
		return *p;
	}

	T& operator()(int x, ...) const
	{
		int p[dim];

		va_list list;
		va_start(list,x);
		p[0] = x;

		for (int i=1; i<dim; i++)
			p[i] = va_arg(list,int);

		va_end(list);

		return operator()(p);
	}


	// buffer I/O
	int buffer_size() const
		{return buffer_size(x0,x1);}

	int buffer_size(const int min[dim], const int max[dim]) const
		{return buffer_size(data,0,min,max);}

	int buffer_size(T* p, int i, const int min[dim], const int max[dim]) const
	{
		int size = 0;
		if (i==dim-1)
			for (int x=min[i]; x<max[i]; x++)
				size += MMSP::buffer_size(*(p+(x-s0[i])*sx[i]));
		else
			for (int x=min[i]; x<max[i]; x++)
				size += buffer_size(p+(x-s0[i])*sx[i],i+1,min,max);
		return size;
	}

	int to_buffer(char* buffer) const
		{return to_buffer(buffer,x0,x1);}

	int to_buffer(char* buffer, const int min[dim], const int max[dim]) const
		{return to_buffer(buffer,data,0,min,max);}

	int to_buffer(char* buffer, T* p, int i, const int min[dim], const int max[dim]) const
	{
		int size = 0;
		if (i==dim-1)
			for (int x=min[i]; x<max[i]; x++)
				size += MMSP::to_buffer(*(p+(x-s0[i])*sx[i]),buffer+size);
		else
			for (int x=min[i]; x<max[i]; x++)
				size += to_buffer(buffer+size,p+(x-s0[i])*sx[i],i+1,min,max);
		return size;
	}

	int from_buffer(char* buffer)
		{return from_buffer(buffer,x0,x1);}

	int from_buffer(char* buffer, const int min[dim], const int max[dim])
		{return from_buffer(buffer,data,0,min,max);}

	int from_buffer(char* buffer, T* p, int i, const int min[dim], const int max[dim])
	{
		int size = 0;
		if (i==dim-1)
			for (int x=min[i]; x<max[i]; x++)
				size += MMSP::from_buffer(*(p+(x-s0[i])*sx[i]),buffer+size);
		else
			for (int x=min[i]; x<max[i]; x++)
				size += from_buffer(buffer+size,p+(x-s0[i])*sx[i],i+1,min,max);
		return size;
	}


	// file I/O
	void input(const char* filename, int GHOSTS=1, int SINGLE=false)
	{
		// file open error check
		std::ifstream input(filename);
		if (!input) {
			std::cerr<<"File input error: could not open ";
			std::cerr<<filename<<"."<<std::endl;
			exit(-1);
		}

		// grid data type error check
		std::string type;
		getline(input,type,'\n');
		if (type!=name(*this)) {
			std::cerr<<"File read error: wrong data type ("<<type<<")."<<std::endl;
			exit(-2);
		}

		// dimension error check
		int dimen;
		input>>dimen;
		if (dimen!=dim) {
			std::cerr<<"File read error: wrong dimension ("<<dimen<<")."<<std::endl;
			exit(-3);
		}

		// read number of fields 
		input>>fields;

		// read grid size
		for (int i=0; i<dim; i++)
			input>>g0[i]>>g1[i];

		// set number of ghosts
		ghosts = 0;

		#ifdef MPI_VERSION
		ghosts = GHOSTS;
		#endif

		// setup grid parameters
		delete [] data;
		setup(SINGLE);

		// read cell spacing
		for (int i=0; i<dim; i++)
			input>>dx[i];

		// ignore trailing endlines
		input.ignore(10,'\n');

		// input grid data
		read(input,GHOSTS);
	}

	void output(const char* filename) const
	{
		#ifndef MPI_VERSION
		// file open error check
		std::ofstream output(filename);
		if (!output) {
			std::cerr<<"File output error: could not open ";
			std::cerr<<filename<<"."<<std::endl;
			exit(-1);
		}

		// write grid data type
		std::string type = name(*this);
		output<<type<<'\n';

		// write grid dimension
		output<<dim<<'\n';

		// write number of fields
		output<<fields<<'\n';

		// write grid size
		for (int i=0; i<dim; i++)
			output<<g0[i]<<" "<<g1[i]<<'\n';

		// write cell spacing
		for (int i=0; i<dim; i++)
			output<<dx[i]<<'\n';
		#endif

		#ifdef MPI_VERSION
		int id = MPI::COMM_WORLD.Get_rank();
		int np = MPI::COMM_WORLD.Get_size();

		std::ofstream output(filename);
		if (id==0) {
			// file open error check
			if (!output) {
				std::cerr<<"File output error: could not open ";
				std::cerr<<filename<<"."<<std::endl;
				exit(-1);
			}

			// write grid data type
			std::string type = name(*this);
			output<<type<<'\n';

			// write grid dimension
			output<<dim<<'\n';

			// write number of fields
			output<<fields<<'\n';

			// write grid size
			for (int i=0; i<dim; i++)
				output<<g0[i]<<" "<<g1[i]<<'\n';

			// write cell spacing
			for (int i=0; i<dim; i++)
				output<<dx[i]<<'\n';
		}
		#endif

		// output grid data
		write(output);
	}

	void read(std::ifstream& file, int GHOSTS=1)
	{
		// read grid data from file
		int lmin[dim], lmax[dim];
		int gmin[dim], gmax[dim];
		for (int i=1; i<dim; i++) {
			lmin[i] = x0[i];
			lmax[i] = x1[i];
			gmin[i] = g0[i];
			gmax[i] = g1[i];
		}

		// iterate through slices in 0th dimension
		for (int x=g0[0]; x<g1[0]; x++) {
			lmin[0] = x;
			lmax[0] = x+1;
			gmin[0] = x;
			gmax[0] = x+1;

			int size;
			file.read(reinterpret_cast<char*>(&size),sizeof(size));
			if (x<x0[0]) file.seekg(size,std::ios::cur);
			else if (x1[0]<=x) break;
			else {
				#ifndef MPI_VERSION
				// read slice data from file
				char* buffer = new char[size];
				file.read(reinterpret_cast<char*>(buffer),size);
				this->from_buffer(buffer,lmin,lmax);
				delete [] buffer;
				#endif

				#ifdef MPI_VERSION
				// read slice data from file
				grid<dim,T> GRID(fields,gmin,gmax,0,true);
				char* buffer = new char[size];
				file.read(reinterpret_cast<char*>(buffer),size);
				GRID.from_buffer(buffer,gmin,gmax);
				delete [] buffer;

				size = GRID.buffer_size(lmin,lmax);
				buffer = new char[size];
				GRID.to_buffer(buffer,lmin,lmax);
				this->from_buffer(buffer,lmin,lmax);
				delete [] buffer;
				#endif
			}
		}
	}

	void write(std::ofstream& file) const
	{
		// write grid data
		int lmin[dim], lmax[dim];
		int gmin[dim], gmax[dim];
		for (int i=1; i<dim; i++) {
			lmin[i] = x0[i];
			lmax[i] = x1[i];
			gmin[i] = g0[i];
			gmax[i] = g1[i];
		}

		// iterate through slices in 0th dimension
		for (int x=g0[0]; x<g1[0]; x++) {
			lmin[0] = x;
			lmax[0] = x+1;
			gmin[0] = x;
			gmax[0] = x+1;

			#ifndef MPI_VERSION
			// write slice data to file
			int size = this->buffer_size(gmin,gmax);
			file.write(reinterpret_cast<const char*>(&size),sizeof(int));
			char* buffer = new char[size];
			this->to_buffer(buffer,gmin,gmax);
			file.write(reinterpret_cast<const char*>(buffer),size);
			delete [] buffer;
			#endif

			#ifdef MPI_VERSION
			int id = MPI::COMM_WORLD.Get_rank();
			int np = MPI::COMM_WORLD.Get_size();

			grid<dim,T> GRID(fields,gmin,gmax,0,true);

			// check if this slice intersects local grid
			int status = 0;
			if (x0[0]<=x and x<x1[0]) status = 1;
			MPI::COMM_WORLD.Send(&status,1,MPI_INT,0,100);

			// if intersection occurs, send data to proc 0
			if (status==1) {
				int size = this->buffer_size(lmin,lmax);
				MPI::COMM_WORLD.Issend(&size,1,MPI_INT,0,200);
				MPI::COMM_WORLD.Issend(lmin,dim,MPI_INT,0,300);
				MPI::COMM_WORLD.Issend(lmax,dim,MPI_INT,0,400);
				char* buffer = new char[size];
				this->to_buffer(buffer,lmin,lmax);
				MPI::COMM_WORLD.Issend(buffer,size,MPI_CHAR,0,500);
				delete [] buffer;
			}

			// prepare global slice in 0th dimension
			if (id==0) {
				// check other processors for data
				for (int i=0; i<np; i++) {
					int status = 0;
					MPI::COMM_WORLD.Recv(&status,1,MPI_INT,i,100);					
					if (status==1) {
						int size;
						int min[dim], max[dim];
						MPI::COMM_WORLD.Recv(&size,1,MPI_INT,i,200);
						MPI::COMM_WORLD.Recv(min,dim,MPI_INT,i,300);
						MPI::COMM_WORLD.Recv(max,dim,MPI_INT,i,400);
						char* buffer = new char[size];
						MPI::COMM_WORLD.Recv(buffer,size,MPI_CHAR,i,500);
						GRID.from_buffer(buffer,min,max);
					}
				}

				// write slice data to file
				int size = GRID.buffer_size();
				file.write(reinterpret_cast<const char*>(&size),sizeof(int));
				char* buffer = new char[size];
				GRID.to_buffer(buffer);
				file.write(reinterpret_cast<const char*>(buffer),size);
				delete [] buffer;
			}

			MPI::COMM_WORLD.Barrier();
			#endif
		}
	}


	// grid parameter "get" functions
	friend int fields(const grid& GRID) {return GRID.fields;}
	friend int ghosts(const grid& GRID) {return GRID.ghosts;}
	friend int x0(const grid& GRID, int i) {return GRID.x0[i];}
	friend int x1(const grid& GRID, int i) {return GRID.x1[i];}
	friend int xmin(const grid& GRID, int i) {return GRID.x0[i];}
	friend int xmax(const grid& GRID, int i) {return GRID.x1[i];}
	friend int xlength(const grid& GRID, int i) {return GRID.x1[i]-GRID.x0[i];}
	friend float dx(const grid& GRID, int i) {return GRID.dx[i];}

	// grid parameters (x direction)
	friend int x0(const grid& GRID) {return GRID.x0[0];}
	friend int x1(const grid& GRID) {return GRID.x1[0];}
	friend int xmin(const grid& GRID) {return GRID.x0[0];}
	friend int xmax(const grid& GRID) {return GRID.x1[0];}
	friend int xlength(const grid& GRID) {return GRID.x1[0]-GRID.x0[0];}
	friend float dx(const grid& GRID) {return GRID.dx[0];}

	// grid parameters (y direction)
	friend int y0(const grid& GRID) {return GRID.x0[1];}
	friend int y1(const grid& GRID) {return GRID.x1[1];}
	friend int ymin(const grid& GRID) {return GRID.x0[1];}
	friend int ymax(const grid& GRID) {return GRID.x1[1];}
	friend int ylength(const grid& GRID) {return GRID.x1[1]-GRID.x0[1];}
	friend float dy(const grid& GRID) {return GRID.dx[1];}

	// grid parameters (z direction)
	friend int z0(const grid& GRID) {return GRID.x0[2];}
	friend int z1(const grid& GRID) {return GRID.x1[2];}
	friend int zmin(const grid& GRID) {return GRID.x0[2];}
	friend int zmax(const grid& GRID) {return GRID.x1[2];}
	friend int zlength(const grid& GRID) {return GRID.x1[2]-GRID.x0[2];}
	friend float dz(const grid& GRID) {return GRID.dx[2];}


	// grid parameter "set" functions
	void set(const char* attribute, int value1[dim], int value2[dim])
	{
		// parallel processes should synchronize data
		#ifdef MPI_VERSION
		int id = MPI::COMM_WORLD.Get_rank();
		int np = MPI::COMM_WORLD.Get_size();

		if (id==0) {
			for (int i=0; i<np; i++) {
				MPI::COMM_WORLD.Send(value1,dim,MPI_INT,i,100);
				MPI::COMM_WORLD.Send(value2,dim,MPI_INT,i,200);
			}
		}

		else {
			MPI::COMM_WORLD.Recv(value1,dim,MPI_INT,0,100);
			MPI::COMM_WORLD.Recv(value2,dim,MPI_INT,0,200);
		}

		MPI::COMM_WORLD.Barrier();
		#endif

		// set attribute
		if (attribute=="boundary") {
			for (int i=0; i<dim; i++) {
				b0[i] = value1[i];
				b1[i] = value1[i];
			}
		}
	}

	void set(const char* attribute, float value1[dim])
	{
		// parallel processes should synchronize data
		#ifdef MPI_VERSION
		int id = MPI::COMM_WORLD.Get_rank();
		int np = MPI::COMM_WORLD.Get_size();

		if (id==0) {
			for (int i=0; i<np; i++) {
				MPI::COMM_WORLD.Send(value1,dim,MPI_FLOAT,i,100);
			}
		}

		else {
			MPI::COMM_WORLD.Recv(value1,dim,MPI_FLOAT,0,100);
		}

		MPI::COMM_WORLD.Barrier();
		#endif

		// set attribute
		if (attribute=="spacing") {
			for (int i=0; i<dim; i++)
				dx[i] = value1[i];
		}
	}

	void set(const char* attribute, int value, ... )
	{
		int value1[dim];
		int value2[dim];

		value1[0] = value;
		value2[0] = value;

		va_list list;
		va_start(list,value);
		for (int i=1; i<dim; i++) {
			value1[i] = va_arg(list,int);
			value2[i] = value1[i];
		}
		va_end(list);

		set(attribute,value1,value2);
	}

	void set(const char* attribute, float value, ... )
	{
		float value1[dim];

		value1[0] = value;

		va_list list;
		va_start(list,value);
		for (int i=1; i<dim; i++) {
			value1[i] = va_arg(list,float);
		}
		va_end(list);

		set(attribute,value1);
	}


	// utility functions
	void swap(grid& GRID)
	{
		// swap grid data 
		T* DATA = data;
		data = GRID.data;
		GRID.data = DATA;

		// swap number of cells
		int CELLS = cells;
		cells = GRID.cells;
		GRID.cells = CELLS;

		// swap number of fields
		int FIELDS = fields;
		fields = GRID.fields;
		GRID.fields = FIELDS;

		// swap number of ghosts 
		int GHOSTS = ghosts;
		ghosts = GRID.ghosts;
		GRID.ghosts = GHOSTS;

		// swap grid parameters
		for (int i=0; i<dim; i++) {
			int G0 = g0[i]; g0[i] = GRID.g0[i]; GRID.g0[i] = G0;
			int G1 = g1[i]; g1[i] = GRID.g1[i]; GRID.g1[i] = G1;
			int X0 = x0[i]; x0[i] = GRID.x0[i]; GRID.x0[i] = X0;
			int X1 = x1[i]; x1[i] = GRID.x1[i]; GRID.x1[i] = X1;
			int S0 = s0[i]; s0[i] = GRID.s0[i]; GRID.s0[i] = S0;
			int S1 = s1[i]; s1[i] = GRID.s1[i]; GRID.s1[i] = S1;
			int SX = sx[i]; sx[i] = GRID.sx[i]; GRID.sx[i] = SX;
			int B0 = b0[i]; b0[i] = GRID.b0[i]; GRID.b0[i] = B0;
			int B1 = b1[i]; b1[i] = GRID.b1[i]; GRID.b1[i] = B1;
			float DX = dx[i]; dx[i] = GRID.dx[i]; GRID.dx[i] = DX;
			int P0 = p0[i]; p0[i] = GRID.p0[i]; GRID.p0[i] = P0;
			int P1 = p1[i]; p1[i] = GRID.p1[i]; GRID.p1[i] = P1;
			int SP = sp[i]; sp[i] = GRID.sp[i]; GRID.sp[i] = SP;
			int N0 = n0[i]; n0[i] = GRID.n0[i]; GRID.n0[i] = N0;
			int N1 = n1[i]; n1[i] = GRID.n1[i]; GRID.n1[i] = N1;
		}
	}

	void copy(const grid& GRID)
	{
		// initialize data
		if (data!=NULL) delete [] data;

		// copy number of cells 
		cells = GRID.cells;

		// copy number of fields
		fields = GRID.fields;

		// copy number of ghosts
		ghosts = GRID.ghosts;

		// copy grid parameters
		for (int i=0; i<dim; i++) {
			g0[i] = GRID.g0[i];
			g1[i] = GRID.g1[i];
			x0[i] = GRID.x0[i];
			x1[i] = GRID.x1[i];
			s0[i] = GRID.s0[i];
			s1[i] = GRID.s1[i];
			sx[i] = GRID.sx[i];
			b0[i] = GRID.b0[i];
			b1[i] = GRID.b1[i];
			dx[i] = GRID.dx[i];
			p0[i] = GRID.p0[i];
			p1[i] = GRID.p1[i];
			sp[i] = GRID.sp[i];
			n0[i] = GRID.n0[i];
			n1[i] = GRID.n1[i];
		}

		// resize data structures 
		data = new T[cells];
		for (int i=0; i<cells; i++)
			resize(data[i],fields);

		// copy grid data 
		for (int i=0; i<cells; i++) {
			int size = MMSP::buffer_size(data[i]);
			char* buffer = new char[size];
			MMSP::to_buffer(GRID.data[i],buffer);
			MMSP::from_buffer(data[i],buffer);
			delete [] buffer;
		}
	}

	void subgrid(grid& GRID, int min[dim], int max[dim]) const
	{
		// copy a subset of grid data to GRID
		int size = this->buffer_size(min,max);
		char* buffer = new char[size];
		this->to_buffer(buffer,min,max);
		GRID.from_buffer(buffer,min,max);
		delete [] buffer;
	}

	void subgrid(grid& GRID, int min0, int max0, ...) const
	{
		int min[dim];
		int max[dim];
		min[0] = min0;
		max[0] = max0;

		va_list list;
		va_start(list,min0);
		for (int i=1; i<dim; i++) {
			min[i] = va_arg(list,int);
			max[i] = va_arg(list,int);
		}
		va_end(list);

		subgrid(GRID,min,max);
	}

	void fill(grid& GRID, int min[dim], int max[dim]) const
	{
		// copy GRID to a subset of grid data
		int size = GRID.buffer_size(min,max);
		char* buffer = new char[size];
		GRID.to_buffer(buffer,min,max);
		this->from_buffer(buffer,min,max);
		delete [] buffer;
	}

	void fill(grid& GRID, int min0, int max0, ...) const
	{
		int min[dim];
		int max[dim];
		min[0] = min0;
		max[0] = max0;

		va_list list;
		va_start(list,min0);
		for (int i=1; i<dim; i++) {
			min[i] = va_arg(list,int);
			max[i] = va_arg(list,int);
		}
		va_end(list);

		fill(GRID,min,max);
	}

	const T& neighbor(int c[dim], int s[dim]) const
	{
		// initialize index
		int p = 0;

		for (int i=0; i<dim; i++) {
			// read arguments
			int x = c[i]+s[i];

			// adjust for boundary conditions
			check_boundary(x,x0[i],x1[i],b0[i],b1[i]);

			// increment position 
			p += (x-s0[i])*sx[i];
		}

		// return cell data
		return data[p];
	}

	const T& neighbor(int c0, int s0, ...) const
	{
		int c[dim];
		int s[dim];
		c[0] = c0;
		s[0] = s0;

		va_list list;
		va_start(list,s0);
		for (int i=1; i<dim; i++) {
			c[i] = va_arg(list,int);
			s[i] = va_arg(list,int);
		}
		va_end(list);

		return neighbor(c,s);
	}


	// parallelization 
	void ghostswap()
	{
		#ifdef MPI_VERSION
		for (int i=0; i<dim; i++) {
			if (1) {
				// send to processor above and receive from processor below
				int send_proc = n1[i];
				int recv_proc = n0[i];

				int send_min[dim], send_max[dim];
				int recv_min[dim], recv_max[dim];
				for (int j=0; j<dim; j++) {
					send_min[j] = x0[j]-ghosts;
					send_max[j] = x1[j]+ghosts;
					recv_min[j] = x0[j]-ghosts;
					recv_max[j] = x1[j]+ghosts;
				}

				send_min[i] = x1[i]-ghosts;
				send_max[i] = x1[i];
				recv_min[i] = x0[i]-ghosts;
				recv_max[i] = x0[i];

				int send_size = this->buffer_size(send_min,send_max);
				int recv_size = 0;
				MPI::COMM_WORLD.Issend(&send_size,1,MPI_INT,send_proc,100);
				MPI::COMM_WORLD.Recv(&recv_size,1,MPI_INT,recv_proc,100);
				char* send_buffer = new char[send_size];
				char* recv_buffer = new char[recv_size];
				this->to_buffer(send_buffer,send_min,send_max);
				MPI::COMM_WORLD.Issend(send_buffer,send_size,MPI_CHAR,send_proc,200);
				MPI::COMM_WORLD.Recv(recv_buffer,recv_size,MPI_CHAR,recv_proc,200);
				MPI::COMM_WORLD.Barrier();
				this->from_buffer(recv_buffer,recv_min,recv_max);
				delete [] send_buffer;
				delete [] recv_buffer;
			}

			if (1) {
				// send to processor below and receive from processor above
				int send_proc = n0[i];
				int recv_proc = n1[i];

				int send_min[dim], send_max[dim];
				int recv_min[dim], recv_max[dim];
				for (int j=0; j<dim; j++) {
					send_min[j] = x0[j]-ghosts;
					send_max[j] = x1[j]+ghosts;
					recv_min[j] = x0[j]-ghosts;
					recv_max[j] = x1[j]+ghosts;
				}

				send_min[i] = x0[i];
				send_max[i] = x0[i]+ghosts;
				recv_min[i] = x1[i];
				recv_max[i] = x1[i]+ghosts;

				int send_size = this->buffer_size(send_min,send_max);
				int recv_size = 0;
				MPI::COMM_WORLD.Issend(&send_size,1,MPI_INT,send_proc,300);
				MPI::COMM_WORLD.Recv(&recv_size,1,MPI_INT,recv_proc,300);
				char* send_buffer = new char[send_size];
				char* recv_buffer = new char[recv_size];
				this->to_buffer(send_buffer,send_min,send_max);
				MPI::COMM_WORLD.Issend(send_buffer,send_size,MPI_CHAR,send_proc,400);
				MPI::COMM_WORLD.Recv(recv_buffer,recv_size,MPI_CHAR,recv_proc,400);
				MPI::COMM_WORLD.Barrier();
				this->from_buffer(recv_buffer,recv_min,recv_max);
				delete [] send_buffer;
				delete [] recv_buffer;
			}
		}
		MPI::COMM_WORLD.Barrier();
		#endif
	}


protected:
	T* data;        // local grid data

	int cells;      // number of cells (including ghosts)
	int fields;     // number of fields
	int ghosts;     // ghost cell depth

	int g0[dim];    // global lower coordinate limit
	int g1[dim];    // global upper coordinate limit

	int x0[dim];    // local lower coordinate limit
	int x1[dim];    // local upper coordinate limit

	int s0[dim];
	int s1[dim];
	int sx[dim];    // local cells in slice 

	int b0[dim];    // boundary condition at x0
	int b1[dim];    // boundary condition at x1

	float dx[dim];  // global cell spacing

	int p0[dim];
	int p1[dim];
	int sp[dim];	// global processors in slice

	int n0[dim];    // neighbor processor at x0
	int n1[dim];    // neighbor processor at x1
};

// buffer I/O functions
template <int dim, typename T> int buffer_size(const grid<dim,T>& GRID)
	{return GRID.buffer_size();}
template <int dim, typename T> int to_buffer(const grid<dim,T>& GRID, char* buffer)
	{return GRID.to_buffer(buffer);}
template <int dim, typename T> int from_buffer(grid<dim,T>& GRID, const char* buffer)
	{return GRID.from_buffer(buffer);}

// file I/O functions
template <int dim, typename T> void read(grid<dim,T>& GRID, std::ifstream& file)
	{GRID.read(file);}
template <int dim, typename T> void write(const grid<dim,T>& GRID, std::ifstream& file)
	{GRID.write(file);}
template <int dim, typename T> void input(grid<dim,T>& GRID, const char* filename, int GHOSTS=1, int SINGLE=false)
	{GRID.input(filename,GHOSTS,SINGLE);}
template <int dim, typename T> void output(const grid<dim,T>& GRID, const char* filename)
	{GRID.output(filename);}

// utility functions
template <int dim, typename T> int length(const grid<dim,T>& GRID) {return 0;}
template <int dim, typename T> void resize(int n, grid<dim,T>& GRID) {}
template <int dim, typename T> void swap(grid<dim,T>& GRID1, grid<dim,T>& GRID2)
	{GRID1.swap(GRID2);}
template <int dim, typename T> void copy(grid<dim,T>& GRID1, grid<dim,T>& GRID2)
	{GRID2.copy(GRID1);}
template <int dim, typename T> std::string name(const grid<dim,T>& GRID) {return std::string("grid:")+name(T());}

template <int dim, typename T> const T& neighbor(grid<dim,T>& GRID, int c[dim], int s[dim])
	{return GRID.neighbor(c,s);}

template <int dim, typename T> const T& neighbor(grid<dim,T>& GRID, int c0, int s0, ...)
	{
		int c[dim];
		int s[dim];
		c[0] = c0;
		s[0] = s0;

		va_list list;
		va_start(list,s0);
		for (int i=1; i<dim; i++) {
			c[i] = va_arg(list,int);
			s[i] = va_arg(list,int);
		}
		va_end(list);

		return GRID.neighbor(c,s);
	}

template <int dim, typename T> void subgrid(const grid<dim,T>& GRID1, grid<dim,T>& GRID2, int min[dim], int max[dim])
	{GRID1.subgrid(GRID2,min,max);}

template <int dim, typename T> void subgrid(const grid<dim,T>& GRID1, grid<dim,T>& GRID2, int min0, int max0, ...)
	{
		int min[dim];
		int max[dim];
		min[0] = min0;
		max[0] = max0;

		va_list list;
		va_start(list,min0);
		for (int i=1; i<dim; i++) {
			min[i] = va_arg(list,int);
			max[i] = va_arg(list,int);
		}
		va_end(list);

		GRID1.subgrid(GRID2,min,max);
	}

template <int dim, typename T> void fill(const grid<dim,T>& GRID1, grid<dim,T>& GRID2, int min[dim], int max[dim])
	{GRID1.fill(GRID2,min,max);}

template <int dim, typename T> void fill(const grid<dim,T>& GRID1, grid<dim,T>& GRID2, int min0, int max0, ...)
	{
		int min[dim];
		int max[dim];
		min[0] = min0;
		max[0] = max0;

		va_list list;
		va_start(list,min0);
		for (int i=1; i<dim; i++) {
			min[i] = va_arg(list,int);
			max[i] = va_arg(list,int);
		}
		va_end(list);

		GRID1.fill(GRID2,min,max);
	}

template <int dim, typename T> void ghostswap(grid<dim,T>& GRID) {GRID.ghostswap();}

template <int dim, typename T> void set(grid<dim,T>& GRID, const char* attribute, int value1[dim], int value2[dim])
	{GRID.set(attribute,value1,value2);}

template <int dim, typename T> void set(grid<dim,T>& GRID, const char* attribute, float value1[dim])
	{GRID.set(attribute,value1);}

template <int dim, typename T> void set(grid<dim,T>& GRID, const char* attribute, int value, ... )
	{
		int value1[dim];
		int value2[dim];

		value1[0] = value;
		value2[0] = value;

		va_list list;
		va_start(list,value);
		for (int i=1; i<dim; i++) {
			value1[i] = va_arg(list,int);
			value2[i] = value1[i];
		}
		va_end(list);

		GRID.set(attribute,value1,value2);
	}

template <int dim, typename T> void set(grid<dim,T>& GRID, const char* attribute, float value, ... )
	{
		int value1[dim];

		value1[0] = value;

		va_list list;
		va_start(list,value);
		for (int i=1; i<dim; i++) {
			value1[i] = va_arg(list,float);
		}
		va_end(list);

		GRID.set(attribute,value1);
	}

template <int dim, typename T> int fields(const grid<dim,T>& GRID) {return fields(GRID);}
template <int dim, typename T> int ghosts(const grid<dim,T>& GRID) {return ghosts(GRID);}
template <int dim, typename T> int x0(const grid<dim,T>& GRID, int i) {return x0(GRID,i);}
template <int dim, typename T> int x1(const grid<dim,T>& GRID, int i) {return x1(GRID,i);}
template <int dim, typename T> int xmin(const grid<dim,T>& GRID, int i) {return xmin(GRID,i);}
template <int dim, typename T> int xmax(const grid<dim,T>& GRID, int i) {return xmax(GRID,i);}
template <int dim, typename T> int xlength(const grid<dim,T>& GRID, int i) {return xlength(GRID,i);}
template <int dim, typename T> float dx(const grid<dim,T>& GRID, int i) {return dx(GRID,i);}

template <int dim, typename T> int x0(const grid<dim,T>& GRID) {return x0(GRID);}
template <int dim, typename T> int x1(const grid<dim,T>& GRID) {return x1(GRID);}
template <int dim, typename T> int xmin(const grid<dim,T>& GRID) {return xmin(GRID);}
template <int dim, typename T> int xmax(const grid<dim,T>& GRID) {return xmax(GRID);}
template <int dim, typename T> int xlength(const grid<dim,T>& GRID) {return xlength(GRID);}
template <int dim, typename T> float dx(const grid<dim,T>& GRID) {return dx(GRID);}

template <int dim, typename T> int y0(const grid<dim,T>& GRID) {return y0(GRID);}
template <int dim, typename T> int y1(const grid<dim,T>& GRID) {return y1(GRID);}
template <int dim, typename T> int ymin(const grid<dim,T>& GRID) {return ymin(GRID);}
template <int dim, typename T> int ymax(const grid<dim,T>& GRID) {return ymax(GRID);}
template <int dim, typename T> int ylength(const grid<dim,T>& GRID) {return ylength(GRID);}
template <int dim, typename T> float dy(const grid<dim,T>& GRID) {return dy(GRID);}

template <int dim, typename T> int z0(const grid<dim,T>& GRID) {return z0(GRID);}
template <int dim, typename T> int z1(const grid<dim,T>& GRID) {return z1(GRID);}
template <int dim, typename T> int zmin(const grid<dim,T>& GRID) {return zmin(GRID);}
template <int dim, typename T> int zmax(const grid<dim,T>& GRID) {return zmax(GRID);}
template <int dim, typename T> int zlength(const grid<dim,T>& GRID) {return zlength(GRID);}
template <int dim, typename T> float dz(const grid<dim,T>& GRID) {return dz(GRID);}

} // namespace MMSP

#endif
