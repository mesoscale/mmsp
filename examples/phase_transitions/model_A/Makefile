# Makefile
# GNU makefile for example "model A" code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -Wall -I $(incdir)
pflags = $(flags) -fpermissive -include mpi.h

# the program
model_A: model_A.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: model_A.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f model_A parallel
