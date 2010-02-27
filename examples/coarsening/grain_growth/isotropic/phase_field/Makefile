# Makefile
# GNU makefile for example grain growth code
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
graingrowth: graingrowth.cpp $(incdir)/PFgrid.hpp PFgrid.isotropy.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: graingrowth.cpp $(incdir)/PFgrid.hpp PFgrid.isotropy.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f graingrowth parallel
