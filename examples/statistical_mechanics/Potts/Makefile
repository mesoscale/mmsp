# Makefile
# GNU makefile for example Potts model code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../include/

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -I $(incdir)
pflags = $(flags) -include mpi.h

# the program
potts: potts.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: potts.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f potts parallel
