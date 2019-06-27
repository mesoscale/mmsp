# hetHC.cpp
# Algorithms for 2D heterogeneous heat conduction, where alpha varies spatially
# Questions/comments to thomas.flint@manchester.ac.uk (Tom Flint)

# includes
incdir = ../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -I $(incdir)
pflags = $(flags) -include mpi.h

# the program
hetHC: hetHC.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: hetHC.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f hetHC parallel
