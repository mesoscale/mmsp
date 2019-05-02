# Makefile
# GNU makefile for grid conversion programs using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# find VTK libraries
VTK_SRC = -I $(shell dirname $(shell locate --limit 1 vtkImageData.h))
VTK_LIB = -L $(shell dirname $(shell locate --limit 1 vtkCommonCore.so))
VTK_LNK = -lvtkCommonCore -lvtkCommonDataModel -lvtkIOCore -lvtkIOXML

# compilers/flags
compiler = g++
flags = -O2 -Wall -I ../include
vtkflags = $(flags) $(VTK_SRC) $(VTK_LIB)
topo = ../algorithms/topology

# conversion programs
utility : mmsp2vti mmsp2pvd mmsp2png mmsp2tsv mmsp2xyz \
          MC2PF MC2sPF \
          PF2MC PF2sPF \
          sPF2MC sPF2PF \
          vox2MC MC2vox \
          mmsp2topo mmsp2z wrongendian

#
#    Utilities for standard visualization formats
#

# convert MMSP grid file to PNG image file type
mmsp2png : mmsp2png.cpp ../include/MMSP.output.h
	$(compiler) $(flags) $< -o $@ -lpng -lz

# convert MMSP grid file to ParaView Data file type
mmsp2pvd : mmsp2pvd.cpp ../include/MMSP.output.h
	$(compiler) $(vtkflags) $< -o $@ $(VTK_LNK) -lz

# convert MMSP grid file to tab-delimited ASCII (TSV) file type
mmsp2tsv : mmsp2tsv.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to compressed VTK Image file type
mmsp2vti : mmsp2vti.cpp ../include/MMSP.output.h
	$(compiler) $(vtkflags) $< -o $@ $(VTK_LNK) -lz

# convert MMSP grid file to XYZ point cloud file type
mmsp2xyz : mmsp2xyz.cpp
	$(compiler) $(flags) $< -o $@ -lz

#
#    Utilities for grid conversions
#

# convert MCgrid file PFgrid file type
MC2PF : MC2PF.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MCgrid file to sPFgrid file type
MC2sPF : MC2sPF.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert PFgrid file MCgrid file type
PF2MC : PF2MC.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert PFgrid file sPFgrid file type
PF2sPF : PF2sPF.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert sPFgrid file MCgrid file type
sPF2MC : sPF2MC.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert sPFgrid file PFgrid file type
sPF2PF : sPF2PF.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MBuilder data to MCgrid file type
vox2MC : vox2MC.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MCgrid file to MBuilder type
MC2vox : MC2vox.cpp
	$(compiler) $(flags) $< -o $@ -lz

#
#    Utilities for niche applications
#

# extract grain boundary topology from MMSP grid file
mmsp2topo : mmsp2topo.cpp
	$(compiler) $(flags) -I $(topo) $< -o $@ -lz

# convert legacy MMSP grid file to compressed MMSP grid file type
mmsp2z : mmsp2z.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file with Little Endian data to Big Endian grid data, and vice versa
# (migrates data between machines with opposite endianness: if you're not sure, you don't need it)
wrongendian : wrongendian.cpp
	$(compiler) $(flags) $< -o $@ -lz -pthread

#
#    Cleanup
#

clean :
	rm -f mmsp2png mmsp2pvd mmsp2tsv mmsp2vti \
	      MC2PF MC2sPF \
	      PF2MC PF2sPF \
	      sPF2MC sPF2PF \
	      vox2MC MC2vox \
	      mmsp2topo mmsp2xyz mmsp2z wrongendian
