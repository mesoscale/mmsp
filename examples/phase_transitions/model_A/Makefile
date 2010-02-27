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
model-A: model-A.cpp model-A.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: model-A.cpp model-A.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f model-A parallel
