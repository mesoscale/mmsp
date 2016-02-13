# If you have a different compiler then you can 
# specify this by override, e.g.: make CC='pgcc'
CC = g++
CCOPTS = -g
PGOPTS = 
LIBS = -lm

core = ../../include

# If you have gnuplot pipes on your system, change this
# value to 'yes' and change the options in the source
# code to get graphical output.
GNUPLOT_I = no

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

ifeq ($(GNUPLOT_I),yes)
	OBJECTS = 1stDiffusion.o gnuplot_i.hpp.gch
else
	OBJECTS = 1stDiffusion.o
endif

#Program dependencies below here.
all : go 
.PHONY : all

go : $(OBJECTS)
	$(CC) $(CCOPTIONS) -o go 1stDiffusion.o $(LIBS)

# Object defs below here.
1stDiffusion.o : 1stDiffusion.cpp
	$(CC) $(CCOPTIONS) -I$(core) -c 1stDiffusion.cpp

gnuplot_i.hpp.gch : gnuplot_i.hpp
	$(CC) $(CCOPTIONS) -c gnuplot_i.hpp


clean :
	rm -rf go *.o *~ *.gch

