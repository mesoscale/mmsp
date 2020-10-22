# Makefile
# GNU makefile for example Cahn-Hilliard model with
# convex splitting and successive over-relaxation
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# includes
incdir = ../../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -Wall -DVANILLA -I $(incdir)
pflags = $(flags) -fpermissive -include mpi.h

# the program
cahn-hilliard: cahn-hilliard.cpp
	$(compiler) $(flags) -fopenmp $< -o $@ -lz

parallel: cahn-hilliard.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f cahn-hilliard parallel
