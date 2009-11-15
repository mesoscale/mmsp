# Makefile
# GNU makefile for example Cahn-Hilliard code
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
cahn-hilliard: cahn-hilliard.cpp cahn-hilliard.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: cahn-hilliard.cpp cahn-hilliard.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f cahn-hilliard parallel
