# Makefile
# GNU makefile for example grid generation
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

# the program
grids: grids.cpp ../include/CAgrid.hpp ../include/MCgrid.hpp ../include/PFgrid.hpp ../include/sparsePF.hpp $(core)
	$(compiler) $(flags) $< -o $@ ; ./grids

clean:
	rm -f grids circle.MC circle.PF circle.sPF sphere.MC sphere.PF sphere.sPF
