# Makefile
# GNU makefile for Poisson solver using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

incdir = ../../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -Wall -I $(incdir)
pflags = $(flags) -fpermissive -include mpi.h

# the program
poisson: poisson.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: poisson.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f poisson parallel
