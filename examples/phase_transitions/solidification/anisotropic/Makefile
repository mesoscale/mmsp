# Makefile
# GNU makefile for dendritic solidification code
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# includes
incdir = ../../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -Wall -I $(incdir)
pflags = $(flags) -fpermissive -include mpi.h

# the program
solidification.out: solidification.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: solidification.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -rf solidification.out parallel
