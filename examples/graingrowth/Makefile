# Makefile
# GNU makefile for example grain growth code
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
graingrowth: graingrowth.cpp ../../include/MCgrid.hpp ../../include/PFgrid.hpp ../../include/sparsePF.hpp $(core) \
             MCgrid.isotropy.hpp PFgrid.isotropy.hpp sparsePF.isotropy.hpp
	$(compiler) $(flags) $< -o $@

parallel: graingrowth.cpp ../../include/MCgrid.hpp ../../include/PFgrid.hpp ../../include/sparsePF.hpp $(core) \
          MCgrid.isotropy.hpp PFgrid.isotropy.hpp sparsePF.isotropy.hpp
	$(pcompiler) $(pflags) $< -o $@

clean:
	rm -f graingrowth parallel
