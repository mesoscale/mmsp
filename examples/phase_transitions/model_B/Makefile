# Makefile
# GNU makefile for example "model B" code
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
model_B: model_B.cpp model_B.hpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

parallel: model_B.cpp model_B.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f model_B parallel
