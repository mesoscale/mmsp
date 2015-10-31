# Makefile
# GNU makefile for example Cahn-Hilliard model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = $(MMSP_PATH)/include

# compilers/flags
compiler = g++
flags = -O3 -Wall -I $(incdir)
gflags = -g -pg -Wall -O2 -I $(incdir)
pcompiler = mpic++
ccompiler = mpic++.mpich2
fcompiler = /usr/mpi/gcc/openmpi-1.8.4rc1/bin/mpic++
pflags = -O3 -Wall -I $(incdir) -include mpi.h
fflags = -O3 -I $(incdir) -include /usr/mpi/gcc/openmpi-1.8.4rc1/include/mpi.h

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp

# the program
cahn-hilliard: cahn-hilliard.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

debug: cahn-hilliard.cpp $(core)
	$(compiler) -DDEBUG $(gflags) $< -o $@ -lz

parallel: cahn-hilliard.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

ctcms: cahn-hilliard.cpp $(core)
	$(ccompiler) $(pflags) $< -o $@ -lz

fnbnd: cahn-hilliard.cpp $(core)
	$(fcompiler) $(fflags) $< -o $@ -lz

clean:
	rm -f cahn-hilliard parallel ctcms debug fnbnd
