# Makefile
# GNU makefile for grid conversion programs using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# compilers/flags
compiler = g++
flags = -O3 -ansi -I ../include/

# dependencies
core = ../include/MMSP.hpp \
       ../include/MMSP.utility.hpp \
       ../include/MMSP.grid.hpp \
       ../include/MMSP.scalar.hpp \
       ../include/MMSP.vector.hpp \
       ../include/MMSP.sparse.hpp 

# conversion programs
utility : mmsp2vtk

mmsp2vtk: mmsp2vtk.cpp $(core)
	$(compiler) $(flags) $< -o $@

clean :
	rm -f mmsp2vtk
