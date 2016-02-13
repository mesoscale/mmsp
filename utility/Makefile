# Makefile
# GNU makefile for grid conversion programs using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# compilers/flags
compiler = g++
flags = -O2 -I ../include/

# dependencies
core = ../include/MMSP.hpp \
       ../include/MMSP.utility.hpp \
       ../include/MMSP.grid.hpp \
       ../include/MMSP.scalar.hpp \
       ../include/MMSP.vector.hpp \
       ../include/MMSP.sparse.hpp 

topo = ../algorithms/topology

# conversion programs
utility : mmsp2vti mmsp2pvd mmsp2png \
          MC2PF MC2sPF \
          PF2MC PF2sPF \
          sPF2MC sPF2PF \
          vox2MC MC2vox \
          mmsp2topo mmsp2z wrongendian

#
#    Utilities for standard visualization formats
#

# convert MMSP grid file to PNG image file type
mmsp2png : mmsp2png.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz -lpng

# convert MMSP grid file to ParaView Data file type
mmsp2pvd : mmsp2pvd.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to VTK Image file type
mmsp2vti : mmsp2vti.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

#
#    Utilities for grid conversions
#

# convert MCgrid file PFgrid file type
MC2PF : MC2PF.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

# convert MCgrid file to sPFgrid file type
MC2sPF : MC2sPF.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

# convert PFgrid file MCgrid file type
PF2MC : PF2MC.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

# convert PFgrid file sPFgrid file type
PF2sPF : PF2sPF.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

# convert sPFgrid file MCgrid file type
sPF2MC : sPF2MC.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

# convert sPFgrid file PFgrid file type
sPF2PF : sPF2PF.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

# convert MBuilder data to MCgrid file type
vox2MC : vox2MC.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

# convert MCgrid file to MBuilder type
MC2vox : MC2vox.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

#
#    Utilities for niche applications
#

# extract grain boundary topology from MMSP grid file
mmsp2topo : mmsp2topo.cpp $(core)
	$(compiler) $(flags) -I $(topo) $< -o $@ -lz

# convert legacy MMSP grid file to compressed MMSP grid file type
mmsp2z : mmsp2z.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file with Little Endian data to Big Endian grid data, and vice versa
# (migrates data between machines with opposite endianness: if you're not sure, you don't need it)
wrongendian : wrongendian.cpp
	$(compiler) $(flags) $< -o $@ -lz -pthread

#
#    Cleanup
#

clean :
	rm -f mmsp2png mmsp2pvd mmsp2vti \
	      MC2PF MC2sPF \
	      PF2MC PF2sPF \
	      sPF2MC sPF2PF \
	      vox2MC MC2vox \
	      mmsp2topo mmsp2z wrongendian
