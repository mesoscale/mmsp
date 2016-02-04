# Makefile
# GNU makefile for single and parallel test programs using MMSP
# and the Google unit test framework (gtest)
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../include
gtdir  = gtest/googletest

# compilers/flags
compiler = g++
flags = -g -I$(incdir)
gflags = -isystem $(gtdir)/include -I$(gtdir) $(flags) -I$(gtdir)/include -pthread
pcompiler = mpic++
pflags = $(flags) -include mpi.h 

# dependencies
core = $(incdir)/MMSP.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp \
       $(incdir)/MMSP.scalar.hpp \
       $(incdir)/MMSP.vector.hpp \
       $(incdir)/MMSP.sparse.hpp 

# the program
test: test.cpp $(core)
	$(compiler) $(flags) $< -o $@

parallel: test.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@

# GTEST


gthdr = $(gtdir)/include/gtest/*.h \
        $(gtdir)/include/gtest/internal/*.h

gtsrc = $(gtdir)/src/*.cc \
        $(gtdir)/src/*.h

gtest-all.o: $(gtsrc) $(gthdr)
	$(compiler) $(gflags) -c $(gtdir)/src/gtest-all.cc

gtest_main.o: $(gtsrc) $(gthdr)
	$(compiler) $(gflags) -c $(gtdir)/src/gtest_main.cc

gtest_main.a: gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

MMSP.unitTestSuite.o: $(core) $(gthdr)
	$(compiler) $(gflags) -c MMSP.unitTestSuite.cpp

unitTestSuite: MMSP.unitTestSuite.o gtest_main.a
	$(compiler) $(gflags) $^ -o $@ -lpthread -lz && \
	./$@

clean:
	rm -f test parallel unitTestSuite *.a *.o
