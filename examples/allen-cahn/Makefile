# Makefile
# GNU makefile for example Allen-Cahn code
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
allen-cahn: allen-cahn.cpp allen-cahn.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: allen-cahn.cpp allen-cahn.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f allen-cahn parallel
