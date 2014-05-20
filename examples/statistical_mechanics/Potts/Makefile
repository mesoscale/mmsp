# Makefile
# GNU makefile for example Potts model code
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
       $(incdir)/MMSP.grid.hpp

# the program
potts: potts.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

parallel: potts.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f potts parallel
