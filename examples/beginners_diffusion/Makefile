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
	OBJECTS = 1stDiffusion.o    gnuplot_i.hpp.gch
	OBJECT1 = MMSPDiffusion.o   gnuplot_i.hpp.gch
	OBJECT2 = MMSPDiffusion2D.o gnuplot_i.hpp.gch
else
	OBJECTS = 1stDiffusion.o
	OBJECT1 = MMSPDiffusion.o
	OBJECT2 = MMSPDiffusion2D.o
endif

#Program dependencies below here.
all : 1stDiffusion MMSPDiffusion MMSPDiffusion2D
.PHONY : all

1stDiffusion : $(OBJECTS)
	$(CC) $(CCOPTIONS) -o $@ 1stDiffusion.o $(LIBS)

MMSPDiffusion : $(OBJECT1)
	$(CC) $(CCOPTIONS) -o $@ MMSPDiffusion.o $(LIBS)

MMSPDiffusion2D : $(OBJECT2)
	$(CC) $(CCOPTIONS) -o $@ MMSPDiffusion2D.o $(LIBS)

# Object defs below here.
1stDiffusion.o : 1stDiffusion.cpp
	$(CC) $(CCOPTIONS) -I$(core) -c $<

MMSPDiffusion.o : MMSPDiffusion.cpp
	$(CC) $(CCOPTIONS) -I$(core) -c $<

MMSPDiffusion2D.o : MMSPDiffusion2D.cpp
	$(CC) $(CCOPTIONS) -I$(core) -c $<

gnuplot_i.hpp.gch : gnuplot_i.hpp
	$(CC) $(CCOPTIONS) -c gnuplot_i.hpp

clean :
	rm -rf 1stDiffusion MMSPDiffusion MMSPDiffusion2D *.o *~ *.gch

