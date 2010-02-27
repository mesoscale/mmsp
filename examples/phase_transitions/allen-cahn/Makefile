# Makefile
# GNU makefile for example Allen-Cahn model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../include

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
allen-cahn: allen-cahn.cpp allen-cahn.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: allen-cahn.cpp allen-cahn.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f allen-cahn parallel
