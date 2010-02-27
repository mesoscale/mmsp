# Makefile
# GNU makefile for example Ostwald ripening code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../../../include

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
ostwald: ostwald.cpp $(incdir)/PFgrid.hpp ostwald.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: ostwald.cpp $(incdir)/PFgrid.hpp ostwald.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f ostwald parallel
