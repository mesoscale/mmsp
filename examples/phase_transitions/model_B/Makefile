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
model-B: model-B.cpp model-B.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: model-B.cpp model-B.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f model-B parallel
