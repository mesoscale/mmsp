# Makefile
# GNU makefile for example Zener pinning code
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../../../../../include

# compilers/flags
compiler = g++
pcompiler = mpic++
flags = -O3 -Wall -I $(incdir)
pflags = $(flags) -fpermissive -include mpi.h

# the program
zener: zener.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: zener.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f zener parallel
