// MMSP.utility.cpp
// Utility function and class implementations for MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include<cstdlib>
#include<iostream>
#include<ciso646>
#include<ctime>

namespace MMSP
{

// MMSP Init function
void Init(int argc, char* argv[])
{
	#ifdef MPI_VERSION
	MPI::Init(argc, argv);
	#endif
}

// MMSP Finalize function
void Finalize()
{
	#ifdef MPI_VERSION
	MPI::Finalize();
	#endif
}

// MMSP Abort function
void Abort(int err)
{
	#ifdef MPI_VERSION
	MPI::COMM_WORLD.Abort(err);
	#endif
	exit(err);
}



// check_boundary: a utility function that adjusts coordinates
// based on limiting coordinates and boundary conditions
void check_boundary(int& x, int x0, int x1, int b0, int b1)
{
	if (x < x0) {
		if (b0 == Neumann or b0 == Dirichlet) x = x0;
		#ifndef MPI_VERSION
		else if (b0 == periodic) x = x1 - (x0 - x);
		#endif
		else if (b0 == mirror) x = 2 * x0 - x;
	} else if (x >= x1) {
		if (b1 == Neumann or b1 == Dirichlet) x = (x1 - 1);
		#ifndef MPI_VERSION
		else if (b1 == periodic) x = x0 + (x - x1);
		#endif
		else if (b1 == mirror) x = 2 * (x1 - 1) - x;
	}
}


// global reducing function
template <typename T> T global(T& value, const char* operation)
{
	// initialize global value
	T global = value;

	#ifdef MPI_VERSION
	int rank = MPI::COMM_WORLD.Get_rank();
	int np = MPI::COMM_WORLD.Get_size();

	if (rank == 0) {
		// receive local values
		for (int i = 1; i < np; i++) {
			T temp;
			int size;
			MPI::COMM_WORLD.Recv(&size, 1, MPI_INT, i, 100);
			char* buffer = new char[size];
			MPI::COMM_WORLD.Recv(buffer, size, MPI_CHAR, i, 200);
			from_buffer(temp, buffer);
			if (buffer != NULL) {
				delete [] buffer;
				buffer=NULL;
			}

			// perform operation
			if (std::string(operation)=="add" or std::string(operation)=="sum")
				global += temp;
			else if (std::string(operation)=="min" or std::string(operation)=="minimum")
				global = min(global, temp);
			else if (std::string(operation)=="max" or std::string(operation)=="maximum")
				global = max(global, temp);
		}

		// send global value
		for (int i = 1; i < np; i++) {
			int size = buffer_size(global);
			MPI::COMM_WORLD.Send(&size, 1, MPI_INT, i, 300);
			char* buffer = new char[size];
			to_buffer(global, buffer);
			MPI::COMM_WORLD.Send(buffer, size, MPI_CHAR, i, 400);
			if (buffer != NULL) {
				delete [] buffer;
				buffer=NULL;
			}
		}
	}

	else {
		// send local value
		int size = buffer_size(value);
		MPI::COMM_WORLD.Send(&size, 1, MPI_INT, 0, 100);
		char* buffer = new char[size];
		to_buffer(value, buffer);
		MPI::COMM_WORLD.Send(buffer, size, MPI_CHAR, 0, 200);
		if (buffer != NULL) {
			delete [] buffer;
			buffer=NULL;
		}

		// receive global value
		MPI::COMM_WORLD.Recv(&size, 1, MPI_INT, 0, 300);
		buffer = new char[size];
		MPI::COMM_WORLD.Recv(buffer, size, MPI_CHAR, 0, 400);
		from_buffer(global, buffer);
		if (buffer != NULL) {
			delete [] buffer;
			buffer=NULL;
		}
	}

	MPI::COMM_WORLD.Barrier();
	#endif

	return global;
}

void print_progress(const int step, const int steps)
{
	/*
	Prints timestamps and a 20-point progress bar to stdout.
	Call once inside the update function (or equivalent).

	for (int step=0; step<steps; step++) {
		if (MPI::COMM_WORLD.Get_rank()==0)
			print_progress(step, steps);
		...
		for (int n=0; n<nodes(grid); n++) {
			...
		}
	}
	*/
	char* timestring;
	static unsigned long tstart;
	struct tm* timeinfo;
	static int iterations = 0;

	if (step==0) {
		iterations++;
		tstart = time(NULL);
		std::time_t rawtime;
		std::time( &rawtime );
		timeinfo = std::localtime( &rawtime );
		timestring = std::asctime(timeinfo);
		timestring[std::strlen(timestring)-1] = '\0';
		std::cout << "No. " << iterations << ":\t" << timestring << " [" << std::flush;
	} else if (step==steps-1) {
		unsigned long deltat = time(NULL)-tstart;
		printf("•] %2luh:%2lum:%2lus",deltat/3600,(deltat%3600)/60,deltat%60);
		std::cout << std::endl;
	} else if ((20 * step) % steps == 0)
		std::cout << "• " << std::flush;
}


} // namespace MMSP
