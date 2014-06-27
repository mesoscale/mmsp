# If you have a different compiler then you can 
# specify this by override, e.g.: make CC='pgcc'
CC = g++
CCOPTS = -g
PGOPTS = 
LIBS = -lm

# If different compilers have different options, then
# those can be handled here.
ifeq ($(CC),gcc)
	CCOPTIONS = $(CCOPTS)
endif

ifeq ($(CC),pgcc)
	CCOPTIONS = $(PGOPTS)
endif

ifeq ($(CC),g++)
	CCOPTIONS = 
endif

#Program dependencies below here.
all : go 
.PHONY : all

go : 1stDiffusion.o gnuplot_i.hpp.gch
	$(CC) $(CCOPTIONS) -o go 1stDiffusion.o $(LIBS)

# Object defs below here.
1stDiffusion.o : 1stDiffusion.cpp
	$(CC) $(CCOPTIONS) -c 1stDiffusion.cpp

gnuplot_i.hpp.gch : gnuplot_i.hpp
	$(CC) $(CCOPTIONS) -c gnuplot_i.hpp

clean :
	rm *.o
	rm *~
	rm *.gch
