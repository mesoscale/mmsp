# Makefile
# GNU makefile for single and parallel test programs using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# includes
incdir = ../include

# compilers/flags
compiler = g++
flags = -O0 -I $(incdir)
pcompiler = mpic++
pflags = -O0 -I $(incdir) -include mpi.h 

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

unitTestSuite: MMSP.unitTestSuite.hpp
	cxxtestgen --error-printer -o MMSP.unitTestSuite.cpp $< && \
	$(compiler) $(flags) --coverage -fprofile-arcs -ftest-coverage -fPIC -o $@ MMSP.unitTestSuite.cpp $(core) && \
	./$@ && \
	gcovr -r . --html --html-details -o unitTestSuite.html

clean:
	rm -f test parallel MMSP.unitTestSuite.cpp unitTestSuite *.html *.gcno *.gcda
