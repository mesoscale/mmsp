# Makefile
# GNU makefile for example "model B" code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -Wall -I $(incdir)
pflags = $(flags) -fpermissive -include mpi.h

# the program
model_B: model_B.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: model_B.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f model_B parallel
