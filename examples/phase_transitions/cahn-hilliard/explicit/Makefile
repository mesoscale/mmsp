# Makefile
# GNU makefile for example Cahn-Hilliard model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -Wall -I $(incdir)
pflags = $(flags) -fpermissive -include mpi.h

# the program
cahn-hilliard: cahn-hilliard.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: cahn-hilliard.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f cahn-hilliard parallel
