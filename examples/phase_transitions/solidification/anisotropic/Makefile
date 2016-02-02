# Makefile
# GNU makefile for dendritic solidification code
# Questions/comments to kellet@rpi.edu (Trevor Keller)

# includes
incdir = ../../../../include

# compilers/flags
compiler = g++ -O3
pcompiler = mpic++ -O3

flags = -I $(incdir)
pflags = $(flags) -include mpi.h

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp \
       $(incdir)/MMSP.sparse.hpp

# the program
solidification.out: solidification.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

mmsp2png: mmsp2png.cpp $(core) /usr/include/IL/devil_cpp_wrapper.hpp
	$(compiler) $(flags) -I /usr/include/IL -include il.h $< -o $@ -lz -lIL -lILU -lILUT

parallel: solidification.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@_solidification.out -lz

clean:
	rm -rf solidification.out parallel_solidification.out mmsp2png
