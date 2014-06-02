# Makefile
# GNU makefile for grid conversion programs using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# compilers/flags
compiler = g++
flags = -O3 -I ../include/ -lz

# dependencies
core = ../include/MMSP.hpp \
       ../include/MMSP.utility.hpp \
       ../include/MMSP.grid.hpp \
       ../include/MMSP.scalar.hpp \
       ../include/MMSP.vector.hpp \
       ../include/MMSP.sparse.hpp 

topo = ../algorithms/topology

# conversion programs
utility : mmsp2vti mmsp2pvd \
          MC2PF MC2sPF \
          PF2MC PF2sPF \
          sPF2MC sPF2PF \
          vox2MC MC2vox \
          TKmmsp2vti mmsp2xyz

mmsp2vti: mmsp2vti.cpp $(core)
	$(compiler) $< -o $@ $(flags)

mmsp2z: mmsp2z.cpp $(core)
	$(compiler) $< -o $@ $(flags)

mmsp2xyz: mmsp2xyz.cpp $(core)
	$(compiler) $< -o $@ $(flags)

mmsp2xyzscalar: mmsp2xyz_scalar.cpp $(core)
	$(compiler) $< -o $@ $(flags)

slicer: mmspslicer.cpp $(core) /usr/include/IL/devil_cpp_wrapper.hpp
	$(compiler) $(flags) -I /usr/include/IL -include il.h $< -o $@ -lz -lIL -lILU -lILUT

TKmmsp2vti: TKmmsp2vti.cpp $(core)
	$(compiler) -I $(topo)/../ $< -o $@ $(flags)

mmsp2topo: mmsp2topo.cpp $(core)
	$(compiler) -I $(topo) $< -o $@ $(flags)

mmsp2pvd: mmsp2pvd.cpp $(core)
	$(compiler) $< -o $@ $(flags)

MC2PF: MC2PF.cpp $(core)
	$(compiler) $< -o $@ $(flags)

MC2sPF: MC2sPF.cpp $(core)
	$(compiler) $< -o $@ $(flags)

PF2MC: PF2MC.cpp $(core)
	$(compiler) $< -o $@ $(flags)

PF2sPF: PF2sPF.cpp $(core)
	$(compiler) $< -o $@ $(flags)

sPF2MC: sPF2MC.cpp $(core)
	$(compiler) $< -o $@ $(flags)

sPF2PF: sPF2PF.cpp $(core)
	$(compiler) $< -o $@ $(flags)

vox2MC: vox2MC.cpp $(core)
	$(compiler) $< -o $@ $(flags)

MC2vox: MC2vox.cpp $(core)
	$(compiler) $< -o $@ $(flags)

clean :
	rm -f TKmmsp2vti slicer mmsp2xyz mmsp2topo mmsp2z mmsp2vti mmsp2pvd MC2PF MC2sPF PF2MC PF2sPF sPF2MC sPF2PF vox2MC MC2vox
