# homHC.cpp
# Algorithms for 2D homogeneous heat conduction
# Questions/comments to thomas.flint@manchester.ac.uk (Tom Flint)

# includes
incdir = ../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -I $(incdir)
pflags = $(flags) -include mpi.h

# the program
homHC: homHC.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: homHC.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f homHC parallel
