# Makefile

# includes
MMSP_PATH=/home/arun/mmsp
incdir = $(MMSP_PATH)/include
utildir = $(MMSP_PATH)/utility
algodir = $(MMSP_PATH)/algorithms

# compilers/flags
compiler = g++ -std=c++11 -O3 
pcompiler = mpic++ -std=c++11 -O3 -Wall -pedantic
flags = -I$(incdir) -I$(algodir) -I$(utildir) 
prefourierflags = -I fftw/api -fopenmp 
postfourierflags = fftw/fftw++.cc  
fourierlinkers = -lfftw3 -lfftw3_omp -lm



# IBM compiler for AMOS
BG_XL = /bgsys/drivers/ppcfloor/comm/xl
BG_INC = -I$(BG_XL)/include
BG_LIB = -L$(BG_XL)/lib
qcompiler = $(BG_XL)/bin/mpixlcxx_r -O3 -qflag=w -qstrict -qmaxmem=-1
qflags = $(BG_INC) $(BG_LIB) $(flags) -I/bgsys/apps/CCNI/zlib/zlib-1.2.7/include -L/bgsys/apps/CCNI/zlib/zlib-1.2.7/lib

# dependencies
core = #$Thermo.hpp \
       #$(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.h \
       $(incdir)/MMSP.grid.h \
       $(incdir)/MMSP.sparse.h \

# the program
sparse: $(core)
	$(compiler) $(flags) $(prefourierflags) main.cpp  $(postfourierflags) $(fourierlinkers) $< -o sparse.out -lz  #-include mpi.h

parallel: $(core)
	$(pcompiler) $(flags) main_parallel.cpp $< -o parallel.out -lz  #-include mpi.h

sparse_serial: $(core)
	$(compiler) $(flags) $(prefourierflags)  main.cpp $(postfourierflags) $(fourierlinkers)  $< -o sparse_serial.out -lz #

bgq: sp-xmpf.cpp $(core)
	$(qcompiler) -DBGQ $(qflags) $< -o q_sparse.out -lz

tool: tool.cpp $(core) /usr/include/IL/devil_cpp_wrapper.hpp
	$(pcompiler) $(flags) -I /usr/include/IL -include il.h $< -o $@ -lz -lIL -lILU -lILUT

mmsp2png : mmsp2png.cpp
	$(compiler) $(flags) $< -o $@ -lz -lpng

# convert MMSP grid file to ParaView Data file type
mmsp2pvd : mmsp2pvd.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to tab-delimited ASCII (TSV) file type
mmsp2tsv : mmsp2tsv.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to VTK Image file type
mmsp2vti : mmsp2vti.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to XYZ point cloud file type
mmsp2xyz : mmsp2xyz.cpp
	$(compiler) $(flags) $< -o $@ -lz

clean:
	rm -rf sparse.out output_* pf_* delta* c_* IE* *~
