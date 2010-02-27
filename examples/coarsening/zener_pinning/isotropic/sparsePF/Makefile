# Makefile
# GNU makefile for example isotropic Zener pinning code
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
       $(incdir)/MMSP.sparse.hpp

# the program
zener: zener.cpp $(incdir)/sparsePF.hpp sparsePF.isotropy.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: zener.cpp $(incdir)/sparsePF.hpp sparsePF.isotropy.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f zener parallel
