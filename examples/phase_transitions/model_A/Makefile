# Makefile
# GNU makefile for example "model A" code
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
model_A: model_A.cpp model_A.hpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

parallel: model_A.cpp model_A.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f model_A parallel
