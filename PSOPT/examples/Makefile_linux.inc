#Usage:  make

SHELL         = /bin/sh
USERHOME      = /home/$(shell whoami)

F2CINC = $(F2CINCLUDE)
F2CLIB = $(F2CLIBRARY)

DMATRIXDIR=../../../dmatrix

CXSPARSE=../../../SuiteSparse/CXSparse
LUSOL=../../../lusol/csrc

prefix = $(USERHOME)/Ipopt-3.12.3
IPOPTINCDIR = -I${prefix}/Ipopt/src/Interfaces -I${prefix}/Ipopt/src/Common -I${prefix}/Ipopt/src/LinAlg
# Directory with libipopt.a
IPOPTLIBDIR = ${exec_prefix}/lib
exec_prefix = ${prefix}

PSOPTDIR    = ../..

PSOPTSRCDIR = ../../src
EXAMPLESDIR = .
LIBDIR        = $(SNOPTDIR)/lib
# SNOPT_LIBS    = $(SNOPTDIR)/lib/libsnopt_cpp.a $(SNOPTDIR)/lib/libsnopt_c.a $(SNOPTDIR)/lib/libsnprint_c.a $(SNOPTDIR)/lib/libblas.a $(SNOPTDIR)/lib/libf2c.a
SNOPT_LIBS    = 
SPARSE_LIBS   = $(LUSOL)/liblusol.a  $(CXSPARSE)/Lib/libcxsparse.a -ldl
# IPOPT_LIBS    =  $(IPOPTLIBDIR)/libipopt.a $(IPOPTLIBDIR)/libcoinhsl.a  $(IPOPTLIBDIR)/ThirdParty/libcoinmumps.a $(IPOPTLIBDIR)/ThirdParty/libcoinmetis.a -lpthread -lgfortran -lblas -llapack -latlas -lf77blas
IPOPT_LIBS = `PKG_CONFIG_PATH=$(prefix)/lib64/pkgconfig:$(prefix)/lib/pkgconfig:$(prefix)/share/pkgconfig: /usr/bin/pkg-config --libs ipopt`
ADOLC_LIBS    = /usr/lib/libadolc.a /usr/lib/libColPack.a
PSOPT_LIBS    = $(PSOPTDIR)/lib/libpsopt.a
DMATRIX_LIBS  = $(DMATRIXDIR)/lib/libdmatrix.a

# LDFLAGS	      = -Wl,--rpath -Wl,$(USERHOME)/Ipopt-3.9.3/lib/coin -L$(IPOPTLIBDIR)  -lm  -ldl    -lgcc_s
LDFLAGS	      = -Wl,--rpath -Wl,$(USERHOME)/Ipopt-3.12.3/lib

CXX           = /usr/bin/g++
CC            = /usr/bin/gcc
CXXFLAGS      = -O0 -g -I$(USERHOME)/adolc_base/include -I$(SNOPTDIR)/cppsrc -I$(DMATRIXDIR)/include -I$(SNOPTDIR)/cppexamples -I$(PSOPTSRCDIR) -DLAPACK -DUNIX -DSPARSE_MATRIX -DUSE_IPOPT -I$(CXSPARSE)/Include -I$(CXSPARSE)/../SuiteSparse_config -I$(LUSOL) $(IPOPTINCDIR) -fomit-frame-pointer -pipe -DNDEBUG -pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -fPIC -DHAVE_MALLOC

CFLAGS        = -O0 -fPIC

# The C++ interface calls some fortran object files, so it is necessary to
# add the Fortran libraries
#
# WARNING: The g95 libraries used here are machine dependent
#
# libg2c and libm are added to satisfy  f2c references

# FLIBS        =   /usr/lib/libf2c.a  -lblas -llapack -lm
FLIBS        =  -lm /usr/lib/libf2c.a  -llapack -lf77blas -lcblas
#FLIBS = ../dmatrix/lib/lapack_LINUX.a ../dmatrix/lib/blas_LINUX.a ../dmatrix/lib/F2CLIBS/libf2c_LINUX.a   -lm

ALL_LIBRARIES = $(PSOPT_LIBS) $(DMATRIX_LIBS) $(FLIBS) $(SPARSE_LIBS) $(IPOPT_LIBS)  $(ADOLC_LIBS)



clean:
	rm -f *.o

$(EXAMPLESDIR)/%.o : $(EXAMPLESDIR)/%.cxx
	$(CXX) -c $(CXXFLAGS) $< -o $@



