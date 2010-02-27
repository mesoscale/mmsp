# Makefile
# GNU makefile for example anisotropic grain growth code
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
       $(incdir)/MMSP.scalar.hpp

# the program
graingrowth: graingrowth.cpp $(incdir)/MCgrid.hpp MCgrid.anisotropy.hpp anisotropy.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: graingrowth.cpp $(incdir)/MCgrid.hpp MCgrid.anisotropy.hpp anisotropy.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f graingrowth parallel
