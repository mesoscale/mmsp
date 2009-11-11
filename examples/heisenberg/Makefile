# Makefile
# GNU makefile for example Heisenberg model code
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
       ../../include/MMSP.scalar.hpp \
       ../../include/MMSP.vector.hpp \
       ../../include/MMSP.sparse.hpp

# the program
heisenberg: heisenberg.cpp heisenberg.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: heisenberg.cpp heisenberg.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f heisenberg parallel
