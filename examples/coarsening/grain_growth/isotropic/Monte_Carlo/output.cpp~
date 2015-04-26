// output.cpp
// Custom output methods for IBM Blue Gene/Q supercomputers:
// specifically the RPI CCI BG/Q, "AMOS".

#ifdef MPI_VERSION

#include<cmath>
#include<sstream>
#include"rdtsc.h"
#include"MMSP.grid.hpp"

namespace MMSP
{

template <int dim,typename T>
double output_bgq(MMSP::grid<dim,T>& GRID, char* filename, int nthreads=1)
{
	/* MPI-IO to the filesystem with writes aligned to blocks */
	MPI::COMM_WORLD.Barrier();
	unsigned long writecycles = rdtsc();
	const unsigned int rank = MPI::COMM_WORLD.Get_rank();
	const unsigned int np = MPI::COMM_WORLD.Get_size();
	MPI_Request request;
	MPI_Status status;
	int mpi_err = 0;
	// Read filesystem block size (using statvfs). Default to 4096 B.
	struct statvfs buf;
	const unsigned long blocksize = (statvfs(".", &buf) == -1)?4096:buf.f_bsize;
	#ifdef DEBUG
	if (rank==0) std::cout<<"Block size is "<<blocksize<<" B."<<std::endl;
	#endif

	// create buffer pointers
	unsigned long* datasizes = NULL;
	unsigned long* offsets = NULL;
	unsigned long* aoffsets = NULL;
	unsigned long* misalignments = NULL;
	char* databuffer=NULL;
	char* headbuffer=NULL;
	char* filebuffer=NULL;
	unsigned int* writeranks=NULL;
	MPI_Request* recvrequests = NULL;
	MPI_Status* recvstatuses = NULL;

	// get grid data to write
	#ifndef RAW
	const unsigned long size=write_buffer_bgq(GRID, databuffer, nthreads);
	#else
	const unsigned long size=write_buffer(GRID, databuffer);
	#endif
	assert(databuffer!=NULL);
	// Generate MMSP header from rank 0
	unsigned long header_offset=0;
	if (rank==0) {
		// get grid data type
		std::string type = name(GRID);

		std::stringstream outstr;
		outstr << type << '\n';
		outstr << dim << '\n';
		outstr << MMSP::fields(GRID) << '\n';

		for (int i=0; i<dim; i++) outstr << MMSP::g0(GRID,i) << " " << MMSP::g1(GRID,i) << '\n'; // global grid dimensions
		for (int i=0; i<dim; i++) outstr << MMSP::dx(GRID,i) << '\n'; // grid spacing

		// Write MMSP header to buffer
		header_offset=outstr.str().size();
		headbuffer = new char[header_offset+sizeof(rank)];
		memcpy(headbuffer, outstr.str().c_str(), header_offset);
		memcpy(headbuffer+header_offset, reinterpret_cast<const char*>(&np), sizeof(np));
		header_offset+=sizeof(rank);
	}
	unsigned long header_offset_send = header_offset;
	MPI::COMM_WORLD.Allreduce(&header_offset_send, &header_offset, 1, MPI_UNSIGNED_LONG, MPI_MAX);
	//	MPI::COMM_WORLD.Bcast(&header_offset, 1, MPI_UNSIGNED_LONG, 0); // broadcast header size from rank 0
	#ifdef DEBUG
	if (rank==0) std::cout<<"Prepared file header."<<std::endl;
	#endif
	MPI::COMM_WORLD.Barrier();
	// Compute file offsets based on buffer sizes
	datasizes = new unsigned long[np];
	MPI::COMM_WORLD.Allgather(&size, 1, MPI_UNSIGNED_LONG, datasizes, 1, MPI_UNSIGNED_LONG);
	#ifdef DEBUG
	if (rank==0) std::cout<<"Synchronized data sizes."<<std::endl;
	#endif
	// Determine disk space requirement
	unsigned long filesize=header_offset;
	for (unsigned int i=0; i<np; ++i) filesize+=datasizes[i];
	MPI::COMM_WORLD.Barrier();

	offsets = new unsigned long[np];
	offsets[0]=header_offset;
	for (unsigned int n=1; n<np; ++n) {
		assert(datasizes[n] < static_cast<unsigned long>(std::numeric_limits<int>::max()));
		offsets[n]=offsets[n-1]+datasizes[n-1];
	}
	offsets[0]=0;
	#ifdef DEBUG
	assert(datasizes[rank]==size);
	if (rank==0) std::cout<<"  Synchronized data offsets on "<<np<<" ranks. Total size: "<<offsets[np-1]+datasizes[np-1]<<" B."<<std::endl;
	#endif
	// Calculate number of  writers & write size
	unsigned long blocks = filesize/blocksize;
	while (blocks*blocksize<filesize)	++blocks;
	const unsigned int nwriters = (blocks>np)?np:blocks;
	const unsigned long writesize=blocksize*(blocks/nwriters);
	assert(writesize % blocksize==0);
	const unsigned long excessblocks=blocks % nwriters;
	bool isWriter=false;
	#ifdef DEBUG
	if (rank==0) std::cout<<"  Preparing "<<nwriters<<" aggregator/writers; writesize is "<<writesize<<" B, with "<<excessblocks<<" excess blocks."<<std::endl;
	#endif
	// Scan to determine which ranks are writers
	writeranks = new unsigned int[nwriters+1];
	aoffsets = new unsigned long[nwriters];
	writeranks[nwriters]=np-1; // generalization for last writer's benefit
	unsigned int temprank=0;
	for (unsigned int w=0; w<nwriters; w++) {
		unsigned long ws=(w<=excessblocks)?writesize+blocksize:writesize;
		// file offset of the w^th writer
		aoffsets[w]=(w>0)?ws+aoffsets[w-1]:0;
		while ((aoffsets[w] > offsets[temprank]+datasizes[temprank]) && temprank<np)
			temprank++;
		writeranks[w]=temprank;
		if (rank==temprank)
			isWriter=true;
		temprank++;
	}
	// Determine which rank to send data to
	unsigned int prevwriter=nwriters, nextwriter=0;
	if (rank==0) {
		prevwriter=0;
	} else {
		while (writeranks[prevwriter]>=rank)
			--prevwriter;
	}
	if (rank>=writeranks[nwriters-1]) {
		nextwriter=nwriters;
	} else {
		while (writeranks[nextwriter]<=rank)
			++nextwriter;
	}
	unsigned long ws = writesize;
	if (nextwriter<=excessblocks)
		ws+=blocksize;
	if (rank>=writeranks[nwriters-1])
		ws=filesize-aoffsets[nwriters-1]; // last block may be only partially-filled

	unsigned long deficiency=0;
	if (rank>0) {
		unsigned long prevws = (prevwriter>=excessblocks)?writesize:writesize+blocksize;
		deficiency = aoffsets[prevwriter]+prevws - offsets[rank];
		if (deficiency>datasizes[rank])
			deficiency=datasizes[rank];
	}
	// Collect block misalignments
	misalignments = new unsigned long[np];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allgather(&deficiency, 1, MPI_UNSIGNED_LONG, misalignments, 1, MPI_UNSIGNED_LONG);
	#ifdef DEBUG
	if (datasizes[rank]-deficiency>ws)
		std::fprintf(stderr, "Error on Rank %u, alignment: buffered %lu B > writesize %lu B.\n", rank, datasizes[rank]-deficiency, ws);
	#endif

	// Accumulate data
	const unsigned int silentranks=writeranks[nextwriter]-rank; // number of MPI ranks between this rank and the next writer
	MPI_Request sendrequest;
	MPI::COMM_WORLD.Barrier();
	if (isWriter) {
		// This rank is a writer.
		assert(misalignments[rank] < datasizes[rank]);
		#ifdef DEBUG
		if (rank>0 && writeranks[prevwriter+1]!=rank)
			std::fprintf(stderr, "Error on Rank %u, writer ID: %u != %u\n", rank, writeranks[prevwriter+1], rank);
		#endif
		// Copy local data into filebuffer
		filebuffer = new char[ws];
		char* p = filebuffer;
		if (rank==0) {
			memcpy(p, headbuffer, header_offset);
			p+=header_offset;
		}
		#ifdef DEBUG
		if (datasizes[rank]-misalignments[rank]>ws)
			std::fprintf(stderr, "Error on Rank %u, memcpy: %lu B > %lu B\n", rank, datasizes[rank]-misalignments[rank], ws);
		#endif
		char* q=databuffer+misalignments[rank];
		memcpy(p, q, datasizes[rank]-misalignments[rank]);
		p+=datasizes[rank]-misalignments[rank];
		// Recv remote data into filebuffer
		if (silentranks>0) {
			recvrequests = new MPI_Request[silentranks];
			recvstatuses = new MPI_Status[silentranks];
		}
		for (unsigned int i=0; i<silentranks && rank+i+1<np; i++) {
			unsigned int recv_proc = rank+i+1;
			assert(recv_proc!=rank && recv_proc<np);
			#ifdef DEBUG
			if (recv_proc<rank || recv_proc>np)
				std::fprintf(stderr, "Error on Rank %u, receiving: recv_proc=%i\n", rank, recv_proc);
			#endif
			unsigned long recv_size = misalignments[recv_proc];
			if (recv_size==0) continue;
			#ifdef DEBUG
			if (p+recv_size>filebuffer+ws)
				std::fprintf(stderr, "Error on Rank %u, receiving from %i: %lu B > %lu B\n", rank, recv_proc, p-filebuffer, ws-recv_size);
			#endif
			MPI_Irecv(p, recv_size, MPI_CHAR, recv_proc, recv_proc, MPI::COMM_WORLD, &recvrequests[i]);
			p+=recv_size;
		}
		#ifdef DEBUG
		if (p-filebuffer!=int(ws))
			std::fprintf(stderr, "Error on Rank %u, total received: %i B != %lu B\n", rank, int(p-filebuffer), ws);
		#endif
		if (rank>0 && misalignments[rank]>0) {
			q=databuffer;
			assert(writeranks[prevwriter]<rank);
			MPI_Isend(q, misalignments[rank], MPI_CHAR, writeranks[prevwriter], rank, MPI::COMM_WORLD, &sendrequest);
		}
	}

	if (misalignments[rank] >= datasizes[rank]) {
		assert(writeranks[prevwriter]<rank && writeranks[prevwriter]<np);
		MPI_Isend(databuffer, datasizes[rank], MPI_CHAR, writeranks[prevwriter], rank, MPI::COMM_WORLD, &sendrequest);
	}
	if (recvrequests != NULL)
		MPI_Waitall(silentranks, recvrequests, recvstatuses);
	if (rank>0) MPI_Wait(&sendrequest, &status);
	MPI::COMM_WORLD.Barrier();

	// file open error check
	#ifdef DEBUG
	if (rank==0) std::cout<<"  Opening "<<std::string(filename)<<" for output."<<std::endl;
	#endif
	MPI_Info info = MPI::INFO_NULL;
	#ifdef BGQ
	MPI_Info_create(&info);
	MPI_Info_set(info, "IBM_largeblock_io", "true");
	#else
	info = MPI::INFO_NULL;
	#endif
	MPI_File output;
	mpi_err = MPI_File_open(MPI::COMM_WORLD, filename, MPI::MODE_WRONLY|MPI::MODE_CREATE, info, &output);
	if (mpi_err != MPI_SUCCESS) {
		char error_string[256];
		int length_of_error_string=256;
		MPI_Error_string(mpi_err, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", rank, error_string);
	}
	if (!output) {
		if (rank==0) std::cerr << "File output error: could not open " << filename << "." << std::endl;
		if (rank==0) std::cerr << "                   If it already exists, delete it and try again." << std::endl;
		exit(-1);
	}
	mpi_err = MPI_File_set_size(output, 0);
	if (mpi_err != MPI_SUCCESS){
		char error_string[256];
		int length_of_error_string=256;
		MPI_Error_string(mpi_err, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", rank, error_string);
	}

	// Write to disk
	if (filebuffer!=NULL) {
		unsigned int w=0;
		while (writeranks[w]!=rank) ++w;
		assert(w<nwriters);
		if (w==nwriters-1)
			assert(filesize-aoffsets[w]==ws);
		mpi_err = MPI_File_iwrite_at(output, aoffsets[w], filebuffer, ws, MPI_CHAR, &request);
		MPI_Wait(&request, &status);
		if (mpi_err != MPI_SUCCESS) {
			char error_string[256];
			int length_of_error_string=256;
			MPI_Error_string(mpi_err, error_string, &length_of_error_string);
			fprintf(stderr, "%3d: %s\n", rank, error_string);
		}
	} else {
		ws = 0; // not a writer
	}

	MPI::COMM_WORLD.Barrier();
	MPI_File_close(&output);
	writecycles = rdtsc() - writecycles;
	MPI_Info_free(&info);
	if (recvrequests!=NULL) {
		delete [] recvrequests;
		recvrequests=NULL;
		delete [] recvstatuses;
		recvstatuses=NULL;
	}
	delete [] misalignments;
	misalignments = NULL;
	delete [] writeranks;
	writeranks=NULL;
	delete [] offsets;
	offsets=NULL;
	delete [] aoffsets;
	aoffsets=NULL;
	delete [] datasizes;
	datasizes=NULL;
	delete [] databuffer;
	databuffer=NULL;
	if (filebuffer!=NULL) {
		delete [] filebuffer;
		filebuffer=NULL;
	}

	return (ws>0)?double(ws)/writecycles:0.0; // bytes per cycle -- needs clock rate info for B/s.
}

template <int dim,typename T>void output_bgq_Pfield(MMSP::grid<dim,T>& GRID, char* filename, int nthreads=1)
{
	/* MPI-IO to the filesystem with writes aligned to blocks */
	MPI::COMM_WORLD.Barrier();
	unsigned long writecycles = rdtsc();
	const unsigned int rank = MPI::COMM_WORLD.Get_rank();
	const unsigned int np = MPI::COMM_WORLD.Get_size();
	MPI_Request request;
	MPI_Status status;
	int mpi_err = 0;
	// Read filesystem block size (using statvfs). Default to 4096 B.
	struct statvfs buf;
	const unsigned long blocksize = (statvfs(".", &buf) == -1)?4096:buf.f_bsize;
	#ifdef DEBUG
	if (rank==0) std::cout<<"Block size is "<<blocksize<<" B."<<std::endl;
	#endif

	// create buffer pointers
	unsigned long* datasizes = NULL;
	unsigned long* offsets = NULL;
	unsigned long* aoffsets = NULL;
	unsigned long* misalignments = NULL;
	char* databuffer=NULL;
	char* headbuffer=NULL;
	char* filebuffer=NULL;
	unsigned int* writeranks=NULL;
	MPI_Request* recvrequests = NULL;
	MPI_Status* recvstatuses = NULL;

	// get grid data to write
	#ifndef RAW
	const unsigned long size=write_buffer_bgq(GRID, databuffer, nthreads);
	#else
if(rank==0) std::cout<<"1111"<<std::endl;
	const unsigned long size=MMSP::write_buffer_Pfield(GRID, databuffer);
if(rank==0) std::cout<<"2222"<<std::endl;
	#endif
	assert(databuffer!=NULL);
	// Generate MMSP header from rank 0
	unsigned long header_offset=0;
	if (rank==0) {
		// get grid data type
		std::string type = name(GRID);

		std::stringstream outstr;
		outstr << type << '\n';
		outstr << dim << '\n';
		outstr << MMSP::fields(GRID) << '\n';

		for (int i=0; i<dim; i++) outstr << MMSP::g0(GRID,i) << " " << MMSP::g1(GRID,i) << '\n'; // global grid dimensions
		for (int i=0; i<dim; i++) outstr << MMSP::dx(GRID,i) << '\n'; // grid spacing

		// Write MMSP header to buffer
		header_offset=outstr.str().size();
		headbuffer = new char[header_offset+sizeof(rank)];
		memcpy(headbuffer, outstr.str().c_str(), header_offset);
		memcpy(headbuffer+header_offset, reinterpret_cast<const char*>(&np), sizeof(np));
		header_offset+=sizeof(rank);
	}
	unsigned long header_offset_send = header_offset;
	MPI::COMM_WORLD.Allreduce(&header_offset_send, &header_offset, 1, MPI_UNSIGNED_LONG, MPI_MAX);
	//	MPI::COMM_WORLD.Bcast(&header_offset, 1, MPI_UNSIGNED_LONG, 0); // broadcast header size from rank 0
	#ifdef DEBUG
	if (rank==0) std::cout<<"Prepared file header."<<std::endl;
	#endif
	MPI::COMM_WORLD.Barrier();
	// Compute file offsets based on buffer sizes
	datasizes = new unsigned long[np];
	MPI::COMM_WORLD.Allgather(&size, 1, MPI_UNSIGNED_LONG, datasizes, 1, MPI_UNSIGNED_LONG);
	#ifdef DEBUG
	if (rank==0) std::cout<<"Synchronized data sizes."<<std::endl;
	#endif
	// Determine disk space requirement
	unsigned long filesize=header_offset;
	for (unsigned int i=0; i<np; ++i) filesize+=datasizes[i];
	MPI::COMM_WORLD.Barrier();

	offsets = new unsigned long[np];
	offsets[0]=header_offset;
	for (unsigned int n=1; n<np; ++n) {
		assert(datasizes[n] < static_cast<unsigned long>(std::numeric_limits<int>::max()));
		offsets[n]=offsets[n-1]+datasizes[n-1];
	}
	offsets[0]=0;
	#ifdef DEBUG
	assert(datasizes[rank]==size);
	if (rank==0) std::cout<<"  Synchronized data offsets on "<<np<<" ranks. Total size: "<<offsets[np-1]+datasizes[np-1]<<" B."<<std::endl;
	#endif
	// Calculate number of  writers & write size
	unsigned long blocks = filesize/blocksize;
	while (blocks*blocksize<filesize)	++blocks;
	const unsigned int nwriters = (blocks>np)?np:blocks;
	const unsigned long writesize=blocksize*(blocks/nwriters);
	assert(writesize % blocksize==0);
	const unsigned long excessblocks=blocks % nwriters;
	bool isWriter=false;
	#ifdef DEBUG
	if (rank==0) std::cout<<"  Preparing "<<nwriters<<" aggregator/writers; writesize is "<<writesize<<" B, with "<<excessblocks<<" excess blocks."<<std::endl;
	#endif
	// Scan to determine which ranks are writers
	writeranks = new unsigned int[nwriters+1];
	aoffsets = new unsigned long[nwriters];
	writeranks[nwriters]=np-1; // generalization for last writer's benefit
	unsigned int temprank=0;
	for (unsigned int w=0; w<nwriters; w++) {
		unsigned long ws=(w<=excessblocks)?writesize+blocksize:writesize;
		// file offset of the w^th writer
		aoffsets[w]=(w>0)?ws+aoffsets[w-1]:0;
		while ((aoffsets[w] > offsets[temprank]+datasizes[temprank]) && temprank<np)
			temprank++;
		writeranks[w]=temprank;
		if (rank==temprank)
			isWriter=true;
		temprank++;
	}
	// Determine which rank to send data to
	unsigned int prevwriter=nwriters, nextwriter=0;
	if (rank==0) {
		prevwriter=0;
	} else {
		while (writeranks[prevwriter]>=rank)
			--prevwriter;
	}
	if (rank>=writeranks[nwriters-1]) {
		nextwriter=nwriters;
	} else {
		while (writeranks[nextwriter]<=rank)
			++nextwriter;
	}
	unsigned long ws = writesize;
	if (nextwriter<=excessblocks)
		ws+=blocksize;
	if (rank>=writeranks[nwriters-1])
		ws=filesize-aoffsets[nwriters-1]; // last block may be only partially-filled

	unsigned long deficiency=0;
	if (rank>0) {
		unsigned long prevws = (prevwriter>=excessblocks)?writesize:writesize+blocksize;
		deficiency = aoffsets[prevwriter]+prevws - offsets[rank];
		if (deficiency>datasizes[rank])
			deficiency=datasizes[rank];
	}
	// Collect block misalignments
	misalignments = new unsigned long[np];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allgather(&deficiency, 1, MPI_UNSIGNED_LONG, misalignments, 1, MPI_UNSIGNED_LONG);
	#ifdef DEBUG
	if (datasizes[rank]-deficiency>ws)
		std::fprintf(stderr, "Error on Rank %u, alignment: buffered %lu B > writesize %lu B.\n", rank, datasizes[rank]-deficiency, ws);
	#endif

	// Accumulate data
	const unsigned int silentranks=writeranks[nextwriter]-rank; // number of MPI ranks between this rank and the next writer
	MPI_Request sendrequest;
	MPI::COMM_WORLD.Barrier();
	if (isWriter) {
		// This rank is a writer.
		assert(misalignments[rank] < datasizes[rank]);
		#ifdef DEBUG
		if (rank>0 && writeranks[prevwriter+1]!=rank)
			std::fprintf(stderr, "Error on Rank %u, writer ID: %u != %u\n", rank, writeranks[prevwriter+1], rank);
		#endif
		// Copy local data into filebuffer
		filebuffer = new char[ws];
		char* p = filebuffer;
		if (rank==0) {
			memcpy(p, headbuffer, header_offset);
			p+=header_offset;
		}
		#ifdef DEBUG
		if (datasizes[rank]-misalignments[rank]>ws)
			std::fprintf(stderr, "Error on Rank %u, memcpy: %lu B > %lu B\n", rank, datasizes[rank]-misalignments[rank], ws);
		#endif
		char* q=databuffer+misalignments[rank];
		memcpy(p, q, datasizes[rank]-misalignments[rank]);
		p+=datasizes[rank]-misalignments[rank];
		// Recv remote data into filebuffer
		if (silentranks>0) {
			recvrequests = new MPI_Request[silentranks];
			recvstatuses = new MPI_Status[silentranks];
		}
		for (unsigned int i=0; i<silentranks && rank+i+1<np; i++) {
			unsigned int recv_proc = rank+i+1;
			assert(recv_proc!=rank && recv_proc<np);
			#ifdef DEBUG
			if (recv_proc<rank || recv_proc>np)
				std::fprintf(stderr, "Error on Rank %u, receiving: recv_proc=%i\n", rank, recv_proc);
			#endif
			unsigned long recv_size = misalignments[recv_proc];
			if (recv_size==0) continue;
			#ifdef DEBUG
			if (p+recv_size>filebuffer+ws)
				std::fprintf(stderr, "Error on Rank %u, receiving from %i: %lu B > %lu B\n", rank, recv_proc, p-filebuffer, ws-recv_size);
			#endif
			MPI_Irecv(p, recv_size, MPI_CHAR, recv_proc, recv_proc, MPI::COMM_WORLD, &recvrequests[i]);
			p+=recv_size;
		}
		#ifdef DEBUG
		if (p-filebuffer!=int(ws))
			std::fprintf(stderr, "Error on Rank %u, total received: %i B != %lu B\n", rank, int(p-filebuffer), ws);
		#endif
		if (rank>0 && misalignments[rank]>0) {
			q=databuffer;
			assert(writeranks[prevwriter]<rank);
			MPI_Isend(q, misalignments[rank], MPI_CHAR, writeranks[prevwriter], rank, MPI::COMM_WORLD, &sendrequest);
		}
	}

	if (misalignments[rank] >= datasizes[rank]) {
		assert(writeranks[prevwriter]<rank && writeranks[prevwriter]<np);
		MPI_Isend(databuffer, datasizes[rank], MPI_CHAR, writeranks[prevwriter], rank, MPI::COMM_WORLD, &sendrequest);
	}
	if (recvrequests != NULL)
		MPI_Waitall(silentranks, recvrequests, recvstatuses);
	if (rank>0) MPI_Wait(&sendrequest, &status);
	MPI::COMM_WORLD.Barrier();

	// file open error check
	#ifdef DEBUG
	if (rank==0) std::cout<<"  Opening "<<std::string(filename)<<" for output."<<std::endl;
	#endif
	MPI_Info info = MPI::INFO_NULL;
	#ifdef BGQ
	MPI_Info_create(&info);
	MPI_Info_set(info, "IBM_largeblock_io", "true");
	#else
	info = MPI::INFO_NULL;
	#endif
	MPI_File output;
	mpi_err = MPI_File_open(MPI::COMM_WORLD, filename, MPI::MODE_WRONLY|MPI::MODE_CREATE, info, &output);
	if (mpi_err != MPI_SUCCESS) {
		char error_string[256];
		int length_of_error_string=256;
		MPI_Error_string(mpi_err, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", rank, error_string);
	}
	if (!output) {
		if (rank==0) std::cerr << "File output error: could not open " << filename << "." << std::endl;
		if (rank==0) std::cerr << "                   If it already exists, delete it and try again." << std::endl;
		exit(-1);
	}
	mpi_err = MPI_File_set_size(output, 0);
	if (mpi_err != MPI_SUCCESS){
		char error_string[256];
		int length_of_error_string=256;
		MPI_Error_string(mpi_err, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", rank, error_string);
	}

	// Write to disk
	if (filebuffer!=NULL) {
		unsigned int w=0;
		while (writeranks[w]!=rank) ++w;
		assert(w<nwriters);
		if (w==nwriters-1)
			assert(filesize-aoffsets[w]==ws);
		mpi_err = MPI_File_iwrite_at(output, aoffsets[w], filebuffer, ws, MPI_CHAR, &request);
		MPI_Wait(&request, &status);
		if (mpi_err != MPI_SUCCESS) {
			char error_string[256];
			int length_of_error_string=256;
			MPI_Error_string(mpi_err, error_string, &length_of_error_string);
			fprintf(stderr, "%3d: %s\n", rank, error_string);
		}
	} else {
		ws = 0; // not a writer
	}

	MPI::COMM_WORLD.Barrier();
	MPI_File_close(&output);
	writecycles = rdtsc() - writecycles;
	MPI_Info_free(&info);
	if (recvrequests!=NULL) {
		delete [] recvrequests;
		recvrequests=NULL;
		delete [] recvstatuses;
		recvstatuses=NULL;
	}
	delete [] misalignments;
	misalignments = NULL;
	delete [] writeranks;
	writeranks=NULL;
	delete [] offsets;
	offsets=NULL;
	delete [] aoffsets;
	aoffsets=NULL;
	delete [] datasizes;
	datasizes=NULL;
	delete [] databuffer;
	databuffer=NULL;
	if (filebuffer!=NULL) {
		delete [] filebuffer;
		filebuffer=NULL;
	}
}

template <int dim,typename T>
double input_bgq(MMSP::grid<dim,T>& GRID, char* filename, int nthreads=1)
{
	/* MPI-IO to the filesystem with writes aligned to blocks */
	MPI::COMM_WORLD.Barrier();
	unsigned long writecycles = rdtsc();
	const unsigned int rank = MPI::COMM_WORLD.Get_rank();
	const unsigned int np = MPI::COMM_WORLD.Get_size();
	MPI_Request request;
	MPI_Status status;
	int mpi_err = 0;

	// Read filesystem block size (using statvfs). Default to 4096 B.
	struct statvfs buf;
	const unsigned long blocksize = (statvfs(".", &buf) == -1)?4096:buf.f_bsize;
	#ifdef DEBUG
	if (rank==0) std::cout<<"Block size is "<<blocksize<<" B."<<std::endl;
	#endif

	// create buffer pointers
	unsigned long* datasizes = NULL;
	unsigned long* offsets = NULL;
	unsigned long* aoffsets = NULL;
	unsigned long* misalignments = NULL;
	char* databuffer=NULL;
	char* headbuffer=NULL;
	char* filebuffer=NULL;
	unsigned int* writeranks=NULL;
	MPI_Request* recvrequests = NULL;
	MPI_Status* recvstatuses = NULL;

	// get grid data to write
	#ifndef RAW
	const unsigned long size=write_buffer_bgq(GRID, databuffer, nthreads);
	#else
	const unsigned long size=write_buffer(GRID, databuffer);
	#endif
	assert(databuffer!=NULL);
	// Generate MMSP header from rank 0
	unsigned long header_offset=0;
	if (rank==0) {
		// get grid data type
		std::string type = name(GRID);

		std::stringstream outstr;
		outstr << type << '\n';
		outstr << dim << '\n';
		outstr << MMSP::fields(GRID) << '\n';

		for (int i=0; i<dim; i++) outstr << MMSP::g0(GRID,i) << " " << MMSP::g1(GRID,i) << '\n'; // global grid dimensions
		for (int i=0; i<dim; i++) outstr << MMSP::dx(GRID,i) << '\n'; // grid spacing

		// Write MMSP header to buffer
		header_offset=outstr.str().size();
		headbuffer = new char[header_offset+sizeof(rank)];
		memcpy(headbuffer, outstr.str().c_str(), header_offset);
		memcpy(headbuffer+header_offset, reinterpret_cast<const char*>(&np), sizeof(np));
		header_offset+=sizeof(rank);
	}
	unsigned long header_offset_send = header_offset;
	MPI::COMM_WORLD.Allreduce(&header_offset_send, &header_offset, 1, MPI_UNSIGNED_LONG, MPI_MAX);
	//	MPI::COMM_WORLD.Bcast(&header_offset, 1, MPI_UNSIGNED_LONG, 0); // broadcast header size from rank 0
	#ifdef DEBUG
	if (rank==0) std::cout<<"Prepared file header."<<std::endl;
	#endif
	MPI::COMM_WORLD.Barrier();

	// Compute file offsets based on buffer sizes
	datasizes = new unsigned long[np];
	MPI::COMM_WORLD.Allgather(&size, 1, MPI_UNSIGNED_LONG, datasizes, 1, MPI_UNSIGNED_LONG);
	#ifdef DEBUG
	if (rank==0) std::cout<<"Synchronized data sizes."<<std::endl;
	#endif

	// Determine disk space requirement
	unsigned long filesize=header_offset;
	for (unsigned int i=0; i<np; ++i) filesize+=datasizes[i];
	MPI::COMM_WORLD.Barrier();

	offsets = new unsigned long[np];
	offsets[0]=header_offset;
	for (unsigned int n=1; n<np; ++n) {
		assert(datasizes[n] < static_cast<unsigned long>(std::numeric_limits<int>::max()));
		offsets[n]=offsets[n-1]+datasizes[n-1];
	}
	offsets[0]=0;
	#ifdef DEBUG
	assert(datasizes[rank]==size);
	if (rank==0) std::cout<<"  Synchronized data offsets on "<<np<<" ranks. Total size: "<<offsets[np-1]+datasizes[np-1]<<" B."<<std::endl;
	#endif

	// Calculate number of  writers & write size
	unsigned long blocks = filesize/blocksize;
	while (blocks*blocksize<filesize)	++blocks;
	const unsigned int nwriters = (blocks>np)?np:blocks;
	const unsigned long writesize=blocksize*(blocks/nwriters);
	assert(writesize % blocksize==0);
	const unsigned long excessblocks=blocks % nwriters;
	bool isWriter=false;
	#ifdef DEBUG
	if (rank==0) std::cout<<"  Preparing "<<nwriters<<" aggregator/writers; writesize is "<<writesize<<" B, with "<<excessblocks<<" excess blocks."<<std::endl;
	#endif

	// Scan to determine which ranks are writers
	writeranks = new unsigned int[nwriters+1];
	aoffsets = new unsigned long[nwriters];
	writeranks[nwriters]=np-1; // generalization for last writer's benefit
	unsigned int temprank=0;
	for (unsigned int w=0; w<nwriters; w++) {
		unsigned long ws=(w<=excessblocks)?writesize+blocksize:writesize;
		// file offset of the w^th writer
		aoffsets[w]=(w>0)?ws+aoffsets[w-1]:0;
		while ((aoffsets[w] > offsets[temprank]+datasizes[temprank]) && temprank<np)
			temprank++;
		writeranks[w]=temprank;
		if (rank==temprank)
			isWriter=true;
		temprank++;
	}

	// Determine which rank to send data to
	unsigned int prevwriter=nwriters, nextwriter=0;
	if (rank==0) {
		prevwriter=0;
	} else {
		while (writeranks[prevwriter]>=rank)
			--prevwriter;
	}
	if (rank>=writeranks[nwriters-1]) {
		nextwriter=nwriters;
	} else {
		while (writeranks[nextwriter]<=rank)
			++nextwriter;
	}

	unsigned long ws = writesize;
	if (nextwriter<=excessblocks)
		ws+=blocksize;
	if (rank>=writeranks[nwriters-1])
		ws=filesize-aoffsets[nwriters-1]; // last block may be only partially-filled

	unsigned long deficiency=0;
	if (rank>0) {
		unsigned long prevws = (prevwriter>=excessblocks)?writesize:writesize+blocksize;
		deficiency = aoffsets[prevwriter]+prevws - offsets[rank];
		if (deficiency>datasizes[rank])
			deficiency=datasizes[rank];
	}
	// Collect block misalignments
	misalignments = new unsigned long[np];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allgather(&deficiency, 1, MPI_UNSIGNED_LONG, misalignments, 1, MPI_UNSIGNED_LONG);

	#ifdef DEBUG
	if (datasizes[rank]-deficiency>ws)
		std::fprintf(stderr, "Error on Rank %u, alignment: buffered %lu B > writesize %lu B.\n", rank, datasizes[rank]-deficiency, ws);
	#endif

	// Accumulate data
	const unsigned int silentranks=writeranks[nextwriter]-rank; // number of MPI ranks between this rank and the next writer
	MPI_Request sendrequest;
	MPI::COMM_WORLD.Barrier();
	if (isWriter) {
		// This rank is a writer.
		assert(misalignments[rank] < datasizes[rank]);
		#ifdef DEBUG
		if (rank>0 && writeranks[prevwriter+1]!=rank)
			std::fprintf(stderr, "Error on Rank %u, writer ID: %u != %u\n", rank, writeranks[prevwriter+1], rank);
		#endif

		// Copy local data into filebuffer
		filebuffer = new char[ws];
		char* p = filebuffer;
		if (rank==0) {
			memcpy(p, headbuffer, header_offset);
			p+=header_offset;
		}
		#ifdef DEBUG
		if (datasizes[rank]-misalignments[rank]>ws)
			std::fprintf(stderr, "Error on Rank %u, memcpy: %lu B > %lu B\n", rank, datasizes[rank]-misalignments[rank], ws);
		#endif
		char* q=databuffer+misalignments[rank];
		memcpy(p, q, datasizes[rank]-misalignments[rank]);
		p+=datasizes[rank]-misalignments[rank];

		// Recv remote data into filebuffer
		if (silentranks>0) {
			recvrequests = new MPI_Request[silentranks];
			recvstatuses = new MPI_Status[silentranks];
		}
		for (unsigned int i=0; i<silentranks && rank+i+1<np; i++) {
			unsigned int recv_proc = rank+i+1;
			assert(recv_proc!=rank && recv_proc<np);
			#ifdef DEBUG
			if (recv_proc<rank || recv_proc>np)
				std::fprintf(stderr, "Error on Rank %u, receiving: recv_proc=%i\n", rank, recv_proc);
			#endif
			unsigned long recv_size = misalignments[recv_proc];
			if (recv_size==0) continue;
			#ifdef DEBUG
			if (p+recv_size>filebuffer+ws)
				std::fprintf(stderr, "Error on Rank %u, receiving from %i: %lu B > %lu B\n", rank, recv_proc, p-filebuffer, ws-recv_size);
			#endif
			MPI_Irecv(p, recv_size, MPI_CHAR, recv_proc, recv_proc, MPI::COMM_WORLD, &recvrequests[i]);
			p+=recv_size;
		}
		#ifdef DEBUG
		if (p-filebuffer!=int(ws))
			std::fprintf(stderr, "Error on Rank %u, total received: %i B != %lu B\n", rank, int(p-filebuffer), ws);
		#endif
		if (rank>0 && misalignments[rank]>0) {
			q=databuffer;
			assert(writeranks[prevwriter]<rank);
			MPI_Isend(q, misalignments[rank], MPI_CHAR, writeranks[prevwriter], rank, MPI::COMM_WORLD, &sendrequest);
		}
	}
	if (misalignments[rank] >= datasizes[rank]) {
		assert(writeranks[prevwriter]<rank && writeranks[prevwriter]<np);
		MPI_Isend(databuffer, datasizes[rank], MPI_CHAR, writeranks[prevwriter], rank, MPI::COMM_WORLD, &sendrequest);
	}
	if (recvrequests != NULL)
		MPI_Waitall(silentranks, recvrequests, recvstatuses);
	if (rank>0) MPI_Wait(&sendrequest, &status);
	MPI::COMM_WORLD.Barrier();

	// file open error check
	#ifdef DEBUG
	if (rank==0) std::cout<<"  Opening "<<std::string(filename)<<" for output."<<std::endl;
	#endif
	MPI_Info info = MPI::INFO_NULL;
	#ifdef BGQ
	MPI_Info_create(&info);
	MPI_Info_set(info, "IBM_largeblock_io", "true");
	#else
	info = MPI::INFO_NULL;
	#endif
	MPI_File output;
	mpi_err = MPI_File_open(MPI::COMM_WORLD, filename, MPI::MODE_WRONLY|MPI::MODE_CREATE, info, &output);
	if (mpi_err != MPI_SUCCESS) {
		char error_string[256];
		int length_of_error_string=256;
		MPI_Error_string(mpi_err, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", rank, error_string);
	}
	if (!output) {
		if (rank==0) std::cerr << "File output error: could not open " << filename << "." << std::endl;
		if (rank==0) std::cerr << "                   If it already exists, delete it and try again." << std::endl;
		exit(-1);
	}
	mpi_err = MPI_File_set_size(output, 0);
	if (mpi_err != MPI_SUCCESS) {
		char error_string[256];
		int length_of_error_string=256;
		MPI_Error_string(mpi_err, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", rank, error_string);
	}

	// Write to disk
	if (filebuffer!=NULL) {
		unsigned int w=0;
		while (writeranks[w]!=rank) ++w;
		assert(w<nwriters);
		if (w==nwriters-1)
			assert(filesize-aoffsets[w]==ws);
		mpi_err = MPI_File_iwrite_at(output, aoffsets[w], filebuffer, ws, MPI_CHAR, &request);
		MPI_Wait(&request, &status);
		if (mpi_err != MPI_SUCCESS) {
			char error_string[256];
			int length_of_error_string=256;
			MPI_Error_string(mpi_err, error_string, &length_of_error_string);
			fprintf(stderr, "%3d: %s\n", rank, error_string);
		}
	} else {
		ws = 0; // not a writer
	}

	MPI::COMM_WORLD.Barrier();
	MPI_File_close(&output);
	writecycles = rdtsc() - writecycles;
	MPI_Info_free(&info);
	if (recvrequests!=NULL) {
		delete [] recvrequests;
		recvrequests=NULL;
		delete [] recvstatuses;
		recvstatuses=NULL;
	}
	delete [] misalignments;
	misalignments = NULL;
	delete [] writeranks;
	writeranks=NULL;
	delete [] offsets;
	offsets=NULL;
	delete [] aoffsets;
	aoffsets=NULL;
	delete [] datasizes;
	datasizes=NULL;
	delete [] databuffer;
	databuffer=NULL;
	if (filebuffer!=NULL) {
		delete [] filebuffer;
		filebuffer=NULL;
	}

	return (ws>0)?double(ws)/writecycles:0.0; // bytes per cycle -- needs clock rate info for B/s.
}

template <int dim,typename T>
double output_split(const MMSP::grid<dim,T>& GRID, char* filename, const int nfiles=1)
{
	// WARNING: THIS FUNCTION HAS NOT BEEN TESTED! IT *WILL* CORRUPT YOUR DATA!

	/* MPI-IO split across multiple files */
	// Function will write a header file named <filename>, plus
	// (nfiles-1) files named <filename%%.dat>.rXXX
	// Assumes that block-alignment DOES NOT MATTER for N>1.


	MPI::COMM_WORLD.Barrier();
	int rank = MPI::COMM_WORLD.Get_rank();
	int np = MPI::COMM_WORLD.Get_size();
	if (nfiles==1) {
		if (rank==0) std::cerr<<"Error in output_split: "<<nfiles<<" specified; should use output_bgq for this case."<<std::endl;
		exit(-1);
	}
	int mpi_err = 0;

	// Set local filename and sub-communicator
	const int ranks_per_file = np/nfiles;
	const int file_number = rank/ranks_per_file;
	int subrank = rank;
	int subnp = np;
	std::string subfile(filename);
	if (file_number>0) {
		// File 0 contains the global header, and uses the fed-in filename
		if (subfile.find_last_of(".") != std::string::npos) {
			subfile = subfile.substr(0,subfile.find_last_of("."))+".";
		}
		std::stringstream nstr;
		nstr << file_number;
		subfile = subfile + "r" + nstr.str();
	}
	char* subfilename = new char[subfile.length()];
	for (int i=0; i<subfile.length(); i++)
		subfilename[i]=subfile[i];

	MPI_Comm subcomm = MPI_COMM_WORLD;
	if (nfiles>1) {
		mpi_err = MPI_Comm_split(MPI::COMM_WORLD, file_number, rank, &subcomm);
		if (mpi_err != MPI_SUCCESS) {
			char error_string[256];
			int length_of_error_string=256;
			MPI_Error_string(mpi_err, error_string, &length_of_error_string);
			fprintf(stderr, "%3d: %s\n", rank, error_string);
		}
		MPI_Comm_rank(subcomm, &subrank);
		MPI_Comm_size(subcomm, &subnp);
	}

	// file open error check
	MPI_Info info = MPI::INFO_NULL;
	#ifdef BGQ
	MPI_Info_create(&info);
	std::string key("IBM_largeblock_io");
	std::string value("true");
	MPI_Info_set(&info, key.c_str(), value.c_str());
	#else
	info = MPI::INFO_NULL;
	#endif
	MPI_File output;
	mpi_err = MPI_File_open(MPI::COMM_WORLD, subfilename, MPI::MODE_WRONLY|MPI::MODE_CREATE, info, &output);
	if (mpi_err != MPI_SUCCESS) {
		char error_string[256];
		int length_of_error_string=256;
		MPI_Error_string(mpi_err, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", rank, error_string);
	}
	if (!output) {
		if (subrank==0) std::cerr << "File output error: could not open " << filename << "." << std::endl;
		if (subrank==0) std::cerr << "                   If it already exists, delete it and try again." << std::endl;
		exit(-1);
	}
	mpi_err = MPI_File_set_size(output, 0);
	if (mpi_err != MPI_SUCCESS) {
		char error_string[256];
		int length_of_error_string=256;
		MPI_Error_string(mpi_err, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", rank, error_string);
	}

	MPI_Request request;
	MPI_Status status;

	unsigned long* databuffer=NULL;

	// Generate global header
	unsigned long header_offset=0;
	char* headbuffer=NULL;
	if (rank==0) {
		// get grid data type
		std::string type = name(GRID);

		std::stringstream outstr;
		outstr << type << '\n';
		outstr << dim << '\n';
		outstr << MMSP::fields(GRID) << '\n';

		for (int i=0; i<dim; i++) outstr << MMSP::g0(GRID,i) << " " << MMSP::g1(GRID,i) << '\n'; // global grid dimensions
		for (int i=0; i<dim; i++) outstr << MMSP::dx(GRID,i) << '\n'; // grid spacing

		// Write MMSP header to file
		header_offset=outstr.str().size();
		headbuffer = new char[header_offset+sizeof(rank)];
		memcpy(headbuffer, outstr.str().c_str(), header_offset);
		memcpy(headbuffer+header_offset, reinterpret_cast<const char*>(&np), sizeof(np));
		header_offset+=sizeof(rank);
	}
	MPI_Barrier(subcomm);
	MPI_Bcast(&header_offset, 1, MPI_UNSIGNED_LONG, 0, subcomm); // broadcast header size from subrank 0; zero in most cases
	#ifdef DEBUG
	if (rank==0) std::cout<<"Prepared file header."<<std::endl;
	#endif
	MPI::COMM_WORLD.Barrier();

	// get grid data to write
	unsigned long size=write_buffer(GRID, databuffer);
	assert(databuffer!=NULL);
	if (rank==0) {
		// Rank 0 holds the global header -- needs to be the first thing written!
		unsigned long* fullbuffer = new unsigned long[size+header_offset];
		memcpy(fullbuffer, headbuffer, header_offset);
		memcpy(fullbuffer+header_offset, databuffer, size);
		delete [] databuffer;
		databuffer = fullbuffer;
		fullbuffer=NULL;
		delete [] headbuffer;
		headbuffer=NULL;
		size+=header_offset;
	}
	assert(databuffer!=NULL);

	// Compute file offsets based on buffer sizes
	unsigned long* datasizes = new unsigned long[subnp];
	MPI_Barrier(subcomm);
	MPI_Allgather(&size, 1, MPI_UNSIGNED_LONG, datasizes, 1, MPI_UNSIGNED_LONG, subcomm);

	// Pre-allocate disk space
	unsigned long filesize=0;
	for (int i=0; i<subnp; ++i) filesize+=datasizes[i];
	MPI_Barrier(subcomm);

	unsigned long* offsets = new unsigned long[subnp];
	offsets[0]=0;
	for (int n=1; n<subnp; ++n) {
		assert(datasizes[n] < static_cast<unsigned long>(std::numeric_limits<int>::max()));
		offsets[n]=offsets[n-1]+datasizes[n-1];
	}
	#ifdef DEBUG
	assert(datasizes[subrank]==size);
	if (subrank==0) std::cout<<"  Synchronized data offsets on "<<subnp<<" ranks. Total size: "<<offsets[subnp-1]+datasizes[subnp-1]<<" B."<<std::endl;
	#endif

	// Write buffer to disk
	MPI_File_sync(output);
	MPI_Barrier(subcomm);
	unsigned long writecycles = 0;
	mpi_err = MPI_File_iwrite_at(output, offsets[subrank], databuffer, datasizes[subrank], MPI_CHAR, &request);
	MPI_Wait(&request, &status);
	writecycles = rdtsc() - writecycles;
	if (mpi_err != MPI_SUCCESS) {
		char error_string[256];
		int length_of_error_string=256;
		MPI_Error_string(mpi_err, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", rank, error_string);
	}
	#ifdef DEBUG
	int error, write_errors=0;
	MPI_Get_count(&status, MPI_INT, &error);
	error++;
	if (error!=1) std::cerr<<"  Error on Rank "<<rank<<": "<<MPI::Get_error_class(error-1)<<std::endl;
	MPI_Allreduce(&error, &write_errors, 1, MPI_INT, MPI_SUM, subcomm);
	assert(write_errors==subnp);
	#endif
	delete [] databuffer;
	databuffer=NULL;

	MPI::COMM_WORLD.Barrier();
	MPI_File_sync(output);
	// Make sure everything's written before closing the file.
	MPI_Offset actual_size;
	MPI_File_get_size(output,&actual_size);
	MPI::COMM_WORLD.Barrier();
	MPI_File_close(&output);

	delete [] offsets;
	offsets=NULL;
	delete [] datasizes;
	datasizes=NULL;

	#ifdef BGQ
	MPI_Info_free(info);
	#endif
	MPI_Comm_free(&subcomm);
	return 0.0;
}


} // namespace MMSP
#endif
