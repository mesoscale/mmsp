# Makefile
# GNU makefile for single and parallel test programs using MMSP
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# includes
incdir = ../include

# compilers/flags
compiler = g++
flags = -g -I$(incdir)
pcompiler = mpic++
pflags = $(flags) -include mpi.h 

# Hello World program
test: test.cpp
	$(compiler) $(flags) $< -o $@

parallel: test.cpp
	$(pcompiler) $(pflags) $< -o $@



# Unit testing depends on the Google C++ testing framework, gtest.
# This is introduced as a git submodule in the gtest sub-directory.
# If yours is a fresh clone of MMSP, you can populate this module
# from the top-level directory by executing
# 	$ git submodule init
# 	$ git submodule update
# Alternatively, clone MMSP recursively:
# 	$ git clone --recursive https://github.com/mesoscale/mmsp.git

gtdir  = gtest/googletest

gflags = -isystem $(gtdir)/include -I$(gtdir) $(flags) \
         -I$(gtdir)/include -pthread

gthdr = $(gtdir)/include/gtest/*.h \
        $(gtdir)/include/gtest/internal/*.h

gtsrc = $(gtdir)/src/*.cc \
        $(gtdir)/src/*.h

# gtest dependencies
gtest-all.o: $(gtsrc) $(gthdr)
	$(compiler) $(gflags) -c $(gtdir)/src/gtest-all.cc

# unit testing program
unitTestSuite: MMSP.unitTestSuite.cpp gtest-all.o
	$(compiler) $(gflags) $^ -o $@ -lpthread -lz

parallelTestSuite: MMSP.unitTestSuite.cpp gtest-all.o
	$(pcompiler) $(gflags) -include mpi.h $^ -o $@ -lpthread -lz

clean:
	rm -f test parallel unitTestSuite parallelTestSuite gtest-all.o
