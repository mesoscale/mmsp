# Makefile
# GNU makefile for example spinodal decomposition code
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
spinodal: spinodal.cpp spinodal.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: spinodal.cpp spinodal.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f spinodal parallel
