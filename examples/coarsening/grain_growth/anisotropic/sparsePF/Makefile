# Makefile
# GNU makefile for example grain growth code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -Wall -I $(incdir)
pflags = $(flags) -fpermissive -include mpi.h

# the program
graingrowth: graingrowth.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: graingrowth.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f graingrowth parallel
