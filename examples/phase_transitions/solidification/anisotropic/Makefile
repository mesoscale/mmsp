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
dendritic.out: dendritic.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

mmsp2png: mmsp2png.cpp $(core) /usr/include/IL/devil_cpp_wrapper.hpp
	$(compiler) $(flags) -I /usr/include/IL -include il.h $< -o $@ -lz -lIL -lILU -lILUT

parallel: dendritic.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@_dendritic.out -lz

clean:
	rm -rf dendritic.out parallel_dendritic.out mmsp2png
