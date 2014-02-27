# Makefile
# GNU makefile for example Heisenberg model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../include/

# compilers/flags
compiler = g++
flags = -O3 -I $(incdir)
pcompiler = mpic++
pflags = -O3 -I $(incdir) -include mpi.h

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp \
       $(incdir)/MMSP.vector.hpp

# the program
heisenberg: heisenberg.cpp heisenberg.hpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

parallel: heisenberg.cpp heisenberg.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f heisenberg parallel
