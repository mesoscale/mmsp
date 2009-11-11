# Makefile
# GNU makefile for example Potts model code
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
potts: potts.cpp potts.hpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: potts.cpp potts.hpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f potts parallel
