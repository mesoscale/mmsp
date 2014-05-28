// wrongendian.cpp
// Change the endianness of MMSP data files
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<zlib.h>
#include<cstring>
#include<cstdlib>
#include<cstdio>
#include<pthread.h>

typedef struct {
	std::ifstream ifile;
	std::ofstream* ofile;
	unsigned long offset;
	int block;
} swap_thread;

template <typename T>
void swap_endian(T& n)
{
	// http://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
	union {
		T u;
		unsigned char u8[sizeof(T)];
	} source, dest;
	source.u = n;
	for (size_t k=0; k<sizeof(T); k++)
		dest.u8[k] = source.u8[sizeof(T)-k-1];
	n = dest.u;
}

template <typename T>
void swap_endian(T* p)
{
	swap_endian<T>(*p);
}

template <typename T>
void swap_buffer(T* buffer, const unsigned long& size)
{
	// Swap endianness of buffer in-place
	T value;
	for (T* p=buffer; p<buffer+size; p+=sizeof(value))
		swap_endian<T>(*p);
}

template <typename T>
void* swap_block_kernel(void* x);

pthread_mutex_t write_lock;

// Define grid variables globablly, to ensure pthreads have access
std::string type;
int dim;
int fields;
int x0[3] = {0, 0, 0};
int x1[3] = {0, 0, 0};
float dx[3] = {1.0, 1.0, 1.0};
int blocks;

bool scalar_type, vector_type, sparse_type;
bool bool_type, char_type, unsigned_char_type, int_type, unsigned_int_type, long_type, unsigned_long_type, short_type, unsigned_short_type, float_type, double_type, long_double_type;

int main(int argc, char* argv[])
{
	// command line error check
	if (argc < 3) {
		std::cout<<"Usage: "<<argv[0]<<" [--help] infile outfile [threads]\n\n";
		exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help" || argc<3) {
		std::cout<<argv[0]<<": Change the endianness of MMSP data files.\n";
		std::cout<<"Usage: "<<argv[0]<<" [--help] infile outfile\n\n";
		std::cout<<"Questions/comments to trevor.keller@gmail.com (Trevor Keller).\n\n";
		exit(0);
	}

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr<<"File input error: could not open "<<argv[1]<<".\n"<<std::endl;
		exit(-1);
	}

	// read data type
	getline(input, type, '\n');

	// grid type error check
	if (type.substr(0, 4) != "grid") {
		std::cerr<<"File input error: file does not contain grid data.\n"<<std::endl;
		exit(-1);
	}

	// parse data type
	scalar_type = (type.find("scalar") != std::string::npos);
	vector_type = (type.find("vector") != std::string::npos);
	sparse_type = (type.find("sparse") != std::string::npos);
	bool_type = (type.find("bool") != std::string::npos);
	char_type = (type.find("char") != std::string::npos);
	unsigned_char_type = (type.find("unsigned char") != std::string::npos);
	int_type = (type.find("int") != std::string::npos);
	unsigned_int_type = (type.find("unsigned int") != std::string::npos);
	long_type = (type.find("long") != std::string::npos);
	unsigned_long_type = (type.find("unsigned long") != std::string::npos);
	short_type = (type.find("short") != std::string::npos);
	unsigned_short_type = (type.find("unsigned short") != std::string::npos);
	float_type = (type.find("float") != std::string::npos);
	double_type = (type.find("double") != std::string::npos);
	long_double_type = (type.find("long double") != std::string::npos);

	if (not bool_type    and
  		not char_type    and  not unsigned_char_type   and
			not int_type     and  not unsigned_int_type    and
			not long_type    and  not unsigned_long_type   and
			not short_type   and  not unsigned_short_type  and
			not float_type   and
			not double_type  and  not long_double_type)
	{
		std::cerr<<"File input error: unknown grid data type ("<<type<<").\n"<<std::endl;
		exit(-1);
	} else if (not sparse_type and not scalar_type and not vector_type) {
		std::cout<<"Warning: "<<type<<" does not contain an MMSP class. Assuming a grid of plain scalars."<<std::endl;
		scalar_type=true;
	}

	// file open error check
	std::ofstream output(argv[2]);
	if (!output) {
		std::cerr<<"File output error: could not open "<<argv[2]<<".\n"<< std::endl;
		exit(-1);
	}

	// Default to 2 threads, unless otherwise specified.
	const unsigned int nthreads=(argc==4)?atoi(argv[3]):2;
	if (nthreads<1) {
		std::cerr<<"POSIX thread error: "<<nthreads<<" threads is too few.\n"<< std::endl;
		exit(-1);
	}

	// write data type
	output<<type<<'\n';

	// copy grid dimension
	input>>dim;
	output<<dim<<'\n';

	// copy number of fields
	input>>fields;
	output<<fields<<'\n';

	// copy grid sizes
	for (int i = 0; i < dim; i++){
		input >> x0[i] >> x1[i];
		output<<x0[i]<<' '<<x1[i]<<'\n';
	}

	// copy cell spacing
	for (int i = 0; i < dim; i++){
		input >> dx[i];
		output<<dx[i]<<'\n';
	}

	// ignore trailing endlines
	input.ignore(10, '\n');


	// copy number of blocks
	input.read(reinterpret_cast<char*>(&blocks), sizeof(blocks));
	swap_endian(blocks);
	output.write(reinterpret_cast<const char*>(&blocks), sizeof(blocks));
	#ifdef DEBUG
	std::cout<<blocks<<" blocks"<<std::endl;
	#endif


	unsigned long pos=input.tellg();
	int b=0;
	while (b < blocks) {
		pthread_t* p_threads = new pthread_t[nthreads];
		swap_thread* swap_threads = new swap_thread[nthreads];
		pthread_attr_t attr;
		pthread_attr_init (&attr);
		int spawned=0;
		for (unsigned int i=0; i<nthreads; i++) {
			if (b<blocks) {
				input.seekg(pos);
				swap_threads[i].block = b-1;
				swap_threads[i].ifile.open(argv[1]);
				swap_threads[i].offset = pos;
				swap_threads[i].ofile = &output;

				pthread_create(&p_threads[i], &attr, swap_block_kernel<float>, (void *) &swap_threads[i] );
				spawned++;

				b++;

				if (b<blocks) {
					pos+=4*dim*sizeof(int)+sizeof(unsigned long);
					unsigned long datasize;
					input.seekg(pos);
					input.read(reinterpret_cast<char*>(&datasize), sizeof(unsigned long));
					swap_endian(datasize);
					pos+=sizeof(unsigned long) + datasize*sizeof(Bytef);
				}
			}
		}
		// Wait for pthreads to exit
		for (int i=0; i<spawned; i++)
			pthread_join(p_threads[i], NULL);
		#ifdef DEBUG
		std::cout<<std::endl;
		#endif

		pthread_attr_destroy(&attr);
		delete [] p_threads;
		delete [] swap_threads;
	}
	#ifdef DEBUG
	std::cout<<"Finished loop."<<std::endl;
	#endif

	input.close();
	output.close();

	std::cout<<"Endianness of "<<argv[1]<<" successfully inverted."<<std::endl;

}

template <typename T>
void* swap_block_kernel(void* x)
{
	swap_thread* st = static_cast<swap_thread*>(x);
	st->ifile.seekg(st->offset);
	// copy block limits
	int lmin[3] = {0, 0, 0};
	int lmax[3] = {0, 0, 0};
	for (int j = 0; j < dim; j++) {
		st->ifile.read(reinterpret_cast<char*>(&lmin[j]), sizeof(lmin[j]));
		st->ifile.read(reinterpret_cast<char*>(&lmax[j]), sizeof(lmax[j]));
		swap_endian<int>(lmin[j]);
		swap_endian<int>(lmax[j]);
	}

	// copy boundary conditions
	int blo[dim];
	int bhi[dim];
	for (int j = 0; j < dim; j++) {
		st->ifile.read(reinterpret_cast<char*>(&blo[j]), sizeof(blo[j]));
		st->ifile.read(reinterpret_cast<char*>(&bhi[j]), sizeof(bhi[j]));
		swap_endian<int>(blo[j]);
		swap_endian<int>(bhi[j]);
	}

	// copy grid data
	unsigned long size_in_mem, size_on_disk;
	st->ifile.read(reinterpret_cast<char*>(&size_in_mem), sizeof(size_in_mem)); // read raw size
	st->ifile.read(reinterpret_cast<char*>(&size_on_disk), sizeof(size_on_disk)); // read compressed size
	swap_endian<unsigned long>(size_in_mem);
	swap_endian<unsigned long>(size_on_disk);
	#ifdef DEBUG
	printf("Block %d: Reading %lu B (%lu KB) into %lu B (%lu KB) buffer.\n", st->block+1, size_on_disk, size_on_disk/1024, size_in_mem, size_in_mem/1024);
	#endif

	// Decompress buffered values
	Bytef* buffer = new Bytef[size_on_disk];
	char* raw=reinterpret_cast<char*>(buffer);
	st->ifile.read(reinterpret_cast<char*>(buffer), size_on_disk);
	st->ifile.close();
	if (size_on_disk!=size_in_mem) {
		size_in_mem=1.125*size_in_mem+12;
		char* raw = new char[size_in_mem];
		// Uncompress data
		int status = uncompress(reinterpret_cast<Bytef*>(raw), &size_in_mem, buffer, size_on_disk);
		if (status<0) std::cerr<<"\nzlib error: "<<status<<std::endl;
		switch( status ) {
			case Z_OK:
				break;
			case Z_MEM_ERROR:
				std::cerr << "Uncompress: out of memory.\n" << std::endl;
				st->ofile->close();
				delete [] raw; raw=NULL;
				if (buffer!=NULL) delete [] buffer; buffer=NULL;
				exit(-1);
				break;
			case Z_BUF_ERROR:
				std::cerr << "Uncompress: output buffer ("<<size_in_mem<<" B) wasn't large enough for data ("<<size_on_disk<<" B).\n" << std::endl;
				st->ofile->close();
				delete [] raw; raw=NULL;
				if (buffer!=NULL) delete [] buffer; buffer=NULL;
				exit(-1);
				break;
		}
		if (buffer!=NULL) delete [] buffer; buffer=NULL;
	}
	// Invert raw data
	if (scalar_type) {
		// MMSP::scalar<T>::to_buffer() copies one T
		char* q=raw;
		if (bool_type or char_type or unsigned_char_type) {
			std::cout<<"Grid type ("<<type<<") does not need byte-swapping."<<std::endl;
			if (raw!=NULL) delete [] raw; raw=NULL;
			if (buffer!=NULL) delete [] buffer; buffer=NULL;
			exit(-1);
		} else if (int_type or unsigned_int_type) {
			while (q<raw+size_in_mem) {
				swap_endian<int>(reinterpret_cast<int*>(q));
				q+=sizeof(int);
			}
		} else if (long_type or unsigned_long_type) {
			while (q<raw+size_in_mem) {
				swap_endian<long>(reinterpret_cast<long*>(q));
				q+=sizeof(long);
			}
		} else if (short_type or unsigned_short_type) {
			while (q<raw+size_in_mem) {
				swap_endian<short>(reinterpret_cast<short*>(q));
				q+=sizeof(short);
			}
		} else if (float_type) {
			while (q<raw+size_in_mem) {
				swap_endian<float>(reinterpret_cast<float*>(q));
				q+=sizeof(float);
			}
		} else if (double_type) {
			while (q<raw+size_in_mem) {
				swap_endian<double>(reinterpret_cast<double*>(q));
				q+=sizeof(double);
			}
		} else if (long_double_type) {
			while (q<raw+size_in_mem) {
				swap_endian<long double>(reinterpret_cast<long double*>(q));
				q+=sizeof(long double);
			}
		} else {
			std::cerr<<"ERROR: Grid type ("<<type<<") is not implemented.\n"<<std::endl;
			if (raw!=NULL) delete [] raw; raw=NULL;
			if (buffer!=NULL) delete [] buffer; buffer=NULL;
			exit(-1);
		}
	} else if (vector_type) {
		// MMSP::vector<T>::to_buffer() copies n * T
		char* q=raw;
		if (bool_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				// sizeof(bool) <= 1B: nothing to swap
				q+=n*sizeof(bool);
			}
		} else if (char_type or unsigned_char_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				// sizeof(char) = 1B: nothing to swap
				q+=n*sizeof(char);
			}
		} else if (int_type or unsigned_int_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				swap_buffer<int>(reinterpret_cast<int*>(q), n);
				q+=n*sizeof(int);
			}
		} else if (long_type or unsigned_long_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				swap_buffer<long>(reinterpret_cast<long*>(q), n);
				q+=n*sizeof(long);
			}
		} else if (short_type or unsigned_short_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				swap_buffer<short>(reinterpret_cast<short*>(q), n);
				q+=n*sizeof(short);
			}
		} else if (float_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				swap_buffer<float>(reinterpret_cast<float*>(q), n);
				q+=n*sizeof(float);
			}
		} else if (double_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				swap_buffer<double>(reinterpret_cast<double*>(q), n);
				q+=n*sizeof(double);
			}
		} else if (long_double_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				swap_buffer<long double>(reinterpret_cast<long double*>(q), n);
				q+=n*sizeof(long double);
			}
		} else {
			std::cerr<<"ERROR: Grid type ("<<type<<") is not implemented.\n"<<std::endl;
			if (raw!=NULL) delete [] raw; raw=NULL;
			if (buffer!=NULL) delete [] buffer; buffer=NULL;
			exit(-1);
		}
	} else if (sparse_type) {
		// MMSP::sparse<T>::to_buffer() copies n * item<T>;
		// each "item" object contains one int and one T
		char* q = raw;
		if (bool_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				for (int j=0; j<n; j++) {
					swap_endian<int>(reinterpret_cast<int*>(q));
					q+=sizeof(int);
					// sizeof(bool) <= 1B: nothing to swap
					q+=sizeof(bool);
				}
			}
		} else if (char_type or unsigned_char_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				for (int j=0; j<n; j++) {
					swap_endian<int>(reinterpret_cast<int*>(q));
					q+=sizeof(int);
					// sizeof(char) = 1B: nothing to swap
					q+=sizeof(char);
				}
			}
		} else if (int_type or unsigned_int_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				for (int j=0; j<n; j++) {
					swap_endian<int>(reinterpret_cast<int*>(q));
					q+=sizeof(int);
					swap_endian<int>(reinterpret_cast<int*>(q));
					q+=sizeof(int);
				}
			}
		} else if (long_type or unsigned_long_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				for (int j=0; j<n; j++) {
					swap_endian<int>(reinterpret_cast<int*>(q));
					q+=sizeof(int);
					swap_endian<long>(reinterpret_cast<long*>(q));
					q+=sizeof(long);
				}
			}
		} else if (short_type or unsigned_short_type) {
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				for (int j=0; j<n; j++) {
					swap_endian<int>(reinterpret_cast<int*>(q));
					q+=sizeof(int);
					swap_endian<short>(reinterpret_cast<short*>(q));
					q+=sizeof(short);
				}
			}
		}	else if (float_type) {
			char* q=raw;
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				for (int j=0; j<n; j++) {
					swap_endian<int>(reinterpret_cast<int*>(q));
					q+=sizeof(int);
					swap_endian<float>(reinterpret_cast<float*>(q));
					q+=sizeof(float);
				}
			}
		}	else if (double_type) {
			char* q=raw;
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				for (int j=0; j<n; j++) {
					swap_endian<int>(reinterpret_cast<int*>(q));
					q+=sizeof(int);
					swap_endian<double>(reinterpret_cast<double*>(q));
					q+=sizeof(double);
				}
			}
		}	else if (long_double_type) {
			char* q=raw;
			while (q<raw+size_in_mem) {
				int n=0;
				swap_endian<int>(reinterpret_cast<int*>(q));
				memcpy(&n, reinterpret_cast<int*>(q), 1);
				q+=sizeof(int);
				for (int j=0; j<n; j++) {
					swap_endian<int>(reinterpret_cast<int*>(q));
					q+=sizeof(int);
					swap_endian<long double>(reinterpret_cast<long double*>(q));
					q+=sizeof(long double);
				}
			}
		} else {
			std::cerr<<"ERROR: Grid type ("<<type<<") is not implemented.\n"<<std::endl;
			if (raw!=NULL) delete [] raw; raw=NULL;
			if (buffer!=NULL) delete [] buffer; buffer=NULL;
			exit(-1);
		}
	} else {
		std::cerr<<"ERROR: Grid type ("<<type<<") is not implemented.\n"<<std::endl;
		if (raw!=NULL) delete [] raw; raw=NULL;
		if (buffer!=NULL) delete [] buffer; buffer=NULL;
		exit(-1);
	}

	// Compress
	size_on_disk=1.25*size_in_mem+12;
	buffer = new Bytef[size_on_disk];
	int status = compress2(buffer, &size_on_disk, reinterpret_cast<const Bytef*>(raw), size_in_mem, 9);
	switch(status) {
		case Z_OK:
			break;
		case Z_MEM_ERROR:
			std::cerr << "Compress: out of memory.\n" << std::endl;
			if (raw!=NULL) delete [] raw; raw=NULL;
			delete [] buffer; buffer=NULL;
			exit(-1);
			break;
		case Z_BUF_ERROR:
			std::cerr << "Compress: output buffer wasn't large enough.\n" << std::endl;
			if (raw!=NULL) delete [] raw; raw=NULL;
			delete [] buffer; buffer=NULL;
			exit(-1);
			break;
	}
	if (raw!=NULL) delete [] raw; raw=NULL;

	pthread_mutex_lock(&write_lock);
	for (int j = 0; j < dim; j++) {
		st->ofile->write(reinterpret_cast<const char*>(&lmin[j]), sizeof(lmin[j]));
		st->ofile->write(reinterpret_cast<const char*>(&lmax[j]), sizeof(lmax[j]));
	}
	for (int j = 0; j < dim; j++) {
		st->ofile->write(reinterpret_cast<const char*>(&blo[j]), sizeof(blo[j]));
		st->ofile->write(reinterpret_cast<const char*>(&bhi[j]), sizeof(bhi[j]));
	}
	st->ofile->write(reinterpret_cast<const char*>(&size_in_mem), sizeof(size_in_mem));
	st->ofile->write(reinterpret_cast<const char*>(&size_on_disk), sizeof(size_on_disk));
	st->ofile->write(reinterpret_cast<const char*>(buffer), size_on_disk);
	pthread_mutex_unlock(&write_lock);

	if (buffer!=NULL) delete [] buffer; buffer=NULL;

	pthread_exit((void*)0);
	return NULL;
}
