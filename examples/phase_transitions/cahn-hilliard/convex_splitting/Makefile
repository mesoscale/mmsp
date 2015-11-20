# Makefile
# GNU makefile for example Cahn-Hilliard model with
# convex splitting and successive over-relaxation
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# includes
incdir = ../../../../include

# compilers/flags
compiler = icc
#flags = -O3 -Wall -std=c++11 -DDEBUG -DVANILLA -I $(incdir)
flags = -O3 -Wall -std=c++11 -DDEBUG -I $(incdir)
gflags = -pg -O2 -Wall -std=c++11 -DDEBUG -DVANILLA -I $(incdir)
pcompiler = mpic++
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

# for profiling
debug: cahn-hilliard.cpp $(core)
	$(compiler) $(gflags) $< -o $@ -lz

clean:
	rm -f cahn-hilliard parallel debug
