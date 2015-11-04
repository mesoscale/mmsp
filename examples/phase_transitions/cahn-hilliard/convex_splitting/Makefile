# Makefile
# GNU makefile for example Cahn-Hilliard model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../../include

# compilers/flags
compiler = icc
gompiler = icc
flags = -O3 -Wall -I $(incdir)
gflags = -g -pg -Wall -O2 -I $(incdir)
pcompiler = mpic++
ccompiler = mpic++.mpich2
pflags = -O3 -Wall -I $(incdir) -include mpi.h

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp

# the program
cahn-hilliard: cahn-hilliard.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

parallel: cahn-hilliard.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

# for profiling and leak checking
debug: cahn-hilliard.cpp $(core)
	$(gompiler) -DDEBUG $(gflags) $< -o $@ -lz

clean:
	rm -f cahn-hilliard parallel debug
