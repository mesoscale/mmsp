# Makefile
# GNU makefile for example anisotropic Zener pinning code
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
zener: zener.cpp zener.hpp anisotropy.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: zener.cpp zener.hpp anisotropy.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f zener parallel
