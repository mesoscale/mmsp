# Makefile
# GNU makefile for grid conversion programs using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# compilers/flags
compiler = g++
flags = -O2 -Wall -I ../include/

# dependencies
core = ../include/MMSP.hpp \
       ../include/MMSP.utility.hpp \
       ../include/MMSP.grid.hpp \
       ../include/MMSP.scalar.hpp \
       ../include/MMSP.vector.hpp \
       ../include/MMSP.sparse.hpp 

topo = ../algorithms/topology

IL_PATH=/users/tnk10/.conda/envs/mmsp

# conversion programs
utility : mmsp2vti mmsp2pvd \
          MC2PF MC2sPF \
          PF2MC PF2sPF \
          sPF2MC sPF2PF \
          vox2MC MC2vox

mmsp2vti: mmsp2vti.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

mmsp2png: mmsp2png.cpp $(core) $(IL_PATH)/include/IL/devil_cpp_wrapper.hpp
	$(compiler) $(flags) -I $(IL_PATH)/include/IL -include il.h $< -o $@ -L$(IL_PATH)/lib -lz -lIL -lILU -lILUT

mmsp2z: mmsp2z.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

slicer: mmspslicer.cpp $(core) $(IL_PATH)/include/IL/devil_cpp_wrapper.hpp
	$(compiler) $(flags) -I $(IL_PATH)/include/IL -include il.h $< -o $@ -L$(IL_PATH)/lib -lz -lIL -lILU -lILUT

mmsp2topo: mmsp2topo.cpp $(core)
	$(compiler) $(flags) -I $(topo) $< -o $@ -lz

mmsp2pvd: mmsp2pvd.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

MC2PF: MC2PF.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

MC2sPF: MC2sPF.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

PF2MC: PF2MC.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

PF2sPF: PF2sPF.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

sPF2MC: sPF2MC.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

sPF2PF: sPF2PF.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

vox2MC: vox2MC.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

MC2vox: MC2vox.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

wrongendian: wrongendian.cpp
	$(compiler) $(flags) $< -o $@ -lz -pthread

clean :
	rm -f slicer mmsp2xyz mmsp2topo mmsp2z mmsp2vti mmsp2pvd MC2PF MC2sPF PF2MC PF2sPF sPF2MC sPF2PF vox2MC MC2vox wrongendian
