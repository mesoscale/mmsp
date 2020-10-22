# Makefile
# GNU makefile for spinodal decomposition example code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -Wall -I $(incdir)
pflags = $(flags) -fpermissive -include mpi.h

# the program
spinodal: spinodal.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: spinodal.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f spinodal parallel
