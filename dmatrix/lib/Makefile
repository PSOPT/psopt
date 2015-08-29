#Usage:  make
#add ../lib to LD_LIBRARY_PATH

SHELL         = /bin/sh



F2CINC = $(F2CINCLUDE)
F2CLIB = $(F2CLIBRARY)

DMATRIXDIR=../../dmatrix
CPPEXAMPLESDIR = $(SNOPTDIR)/cppexamples
CXSPARSE=../../SuiteSparse/CXSparse
LUSOL=../../lusol/csrc


# Directory with header files
IPOPTINCDIR = ${prefix}/include/coin
# Directory with libipopt.a
IPOPTLIBDIR = ${exec_prefix}/lib
exec_prefix = ${prefix}


PSOPTSRCDIR = ../src


CXX           = /usr/bin/g++
CC            = /usr/bin/gcc
CXXFLAGS      = -O0 -g -Wno-long-long -I$(SNOPTDIR)/cppsrc -I$(DMATRIXDIR)/include -I$(SNOPTDIR)/cppexamples -I$(PSOPTSRCDIR) -DLAPACK -DUNIX -DSPARSE_MATRIX -DUSE_SNOPT -DUSE_IPOPT -I$(CXSPARSE)/Include -I$(LUSOL) -I$(CXSPARSE)/../SuiteSparse_config -I$(IPOPTINCDIR) -fomit-frame-pointer -pipe -DNDEBUG -pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -fPIC -DHAVE_MALLOC

CFLAGS        = -O0 -fPIC


DMATRIXLIB = libdmatrix.a

$(DMATRIXLIB):  $(DMATRIXLIB)($(DMATRIXDIR)/src/dmatrixv.o)


$(DMATRIXDIR)/src/%.o : $(DMATRIXDIR)/src/%.cxx $(DMATRIXDIR)/include/%.h
	$(CXX) -c $(CXXFLAGS) $< -o $@


clean:

	(rm -f *.o;rm -f ../src/*.o)

distclean:

	(rm -f *.a)
