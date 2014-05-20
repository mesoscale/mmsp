# Makefile
# GNU makefile for example Zener pinning code
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
zener: zener.cpp anisotropy.hpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

parallel: zener.cpp anisotropy.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f zener parallel
