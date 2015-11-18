# Makefile
# GNU makefile for example Cahn-Hilliard model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../../include

# compilers/flags
compiler = icc
gompiler = icc
flags = -O3 -Wall -ffast-math -std=c++11 -DDEBUG -DVANILLA -I $(incdir)
#flags = -O3 -Wall -ffast-math -std=c++11 -DDEBUG -I $(incdir)
gflags = -g -pg -Wall -O2 -I $(incdir)
pcompiler = mpic++
ccompiler = mpic++.mpich2
pflags = $(flags) -include mpi.h

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
	$(gompiler) $(gflags) $< -o $@ -lz

tarball:
	tar -cvzf cahn-hilliard.tar.gz cahn-hilliard.cpp cahn-hilliard.hpp Makefile notes

clean:
	rm -f cahn-hilliard parallel debug cahn-hilliard.tar.gz
