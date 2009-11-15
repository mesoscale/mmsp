# Makefile
# GNU makefile for example "model B" code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# compilers/flags
compiler = g++
flags = -O3 -ansi -I ../../include/
pcompiler = mpic++
pflags = -O3 -ansi -I ../../include/ -include mpi.h

# dependencies
core = ../../include/MMSP.hpp \
       ../../include/MMSP.utility.hpp \
       ../../include/MMSP.grid.hpp \
       ../../include/MMSP.scalar.hpp

# the program
model-B: model-B.cpp model-B.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: model-B.cpp model-B.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f model-B parallel
