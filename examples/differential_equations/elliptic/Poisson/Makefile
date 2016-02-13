# Makefile
# GNU makefile for Poisson solver using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

incdir = ../../../../include

# compilers/flags
compiler = g++
flags = -O3 -ansi -I $(incdir)
pcompiler = mpic++
pflags = -O3 -ansi -I $(incdir) -include mpi.h 

# dependencies
core = $(incdir)/MMSP.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp

# the program
poisson: poisson.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

parallel: poisson.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f poisson parallel
