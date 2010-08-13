# Makefile
# GNU makefile for grid conversion programs using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# compilers/flags
compiler = g++
flags = -O3 -I ../include/

# dependencies
core = ../include/MMSP.hpp \
       ../include/MMSP.utility.hpp \
       ../include/MMSP.grid.hpp \
       ../include/MMSP.scalar.hpp \
       ../include/MMSP.vector.hpp \
       ../include/MMSP.sparse.hpp 

# conversion programs
utility : mmsp2vti mmsp2pvd \
          MC2PF MC2sPF \
          PF2MC PF2sPF \
          sPF2MC sPF2PF \
          vox2MC MC2vox

mmsp2vti: mmsp2vti.cpp $(core)
	$(compiler) $(flags) $< -o $@

mmsp2pvd: mmsp2pvd.cpp $(core)
	$(compiler) $(flags) $< -o $@

MC2PF: MC2PF.cpp $(core)
	$(compiler) $(flags) $< -o $@

MC2sPF: MC2sPF.cpp $(core)
	$(compiler) $(flags) $< -o $@

PF2MC: PF2MC.cpp $(core)
	$(compiler) $(flags) $< -o $@

PF2sPF: PF2sPF.cpp $(core)
	$(compiler) $(flags) $< -o $@

sPF2MC: sPF2MC.cpp $(core)
	$(compiler) $(flags) $< -o $@

sPF2PF: sPF2PF.cpp $(core)
	$(compiler) $(flags) $< -o $@

vox2MC: vox2MC.cpp $(core)
	$(compiler) $(flags) $< -o $@

MC2vox: MC2vox.cpp $(core)
	$(compiler) $(flags) $< -o $@

clean :
	rm -f mmsp2vti mmsp2pvd MC2PF MC2sPF PF2MC PF2sPF sPF2MC sPF2PF vox2MC MC2vox
