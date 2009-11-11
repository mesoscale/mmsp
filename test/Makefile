# Makefile
# GNU makefile for single and parallel test programs using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# compilers/flags
compiler = g++
flags = -O3 -ansi -I ../include/
pcompiler = mpic++
pflags = -O3 -ansi -I ../include/ -include mpi.h 

# dependencies
core = ../include/MMSP.hpp \
       ../include/MMSP.utility.hpp \
       ../include/MMSP.grid.hpp \
       ../include/MMSP.scalar.hpp \
       ../include/MMSP.vector.hpp \
       ../include/MMSP.sparse.hpp 

# the program
test: test.cpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: test.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f test parallel
