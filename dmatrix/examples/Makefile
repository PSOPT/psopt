
#Usage:  make all

SHELL         = /bin/sh
F2CINC = $(F2CINCLUDE)
F2CLIB = $(F2CLIBRARY)
DMATRIXDIR= ../../dmatrix
CXSPARSE=../../SuiteSparse/CXSparse
LUSOL=../../lusol/csrc



EXAMPLESDIR = .

SPARSE_LIBS   = $(LUSOL)/liblusol.a  $(CXSPARSE)/Lib/libcxsparse.a 

DMATRIX_LIBS  = $(DMATRIXDIR)/lib/libdmatrix.a 

FLIBS        =   -llapack -lblas /usr/lib/libf2c.a

LDFLAGS	      = -Wl,--rpath  -lgcc_s

ALL_LIBRARIES =   $(DMATRIX_LIBS) $(FLIBS) $(SPARSE_LIBS) -ldl -lm -lgcc_s

CXX           = /usr/bin/g++
CC            = /usr/bin/gcc
CXXFLAGS      = -O0 -g -I$(DMATRIXDIR)/include -DLAPACK -DUNIX -DSPARSE_MATRIX  -I$(CXSPARSE)/Include -I$(CXSPARSE)/../SuiteSparse_config -I$(LUSOL) -fomit-frame-pointer -pipe -DNDEBUG -pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -fPIC -DHAVE_MALLOC

CFLAGS        = -O0 -fPIC


RICC  = ricc   
COLON = colon
EXPM  = expm
RAND  = rand
SPARSE = sparse


RICC_O = $(RICC:%=$(EXAMPLESDIR)/%.o)
COLON_O  = $(COLON:%=$(EXAMPLESDIR)/%.o)
EXPM_O  = $(EXPM:%=$(EXAMPLESDIR)/%.o)
RAND_O  = $(RAND:%=$(EXAMPLESDIR)/%.o)
SPARSE_O  = $(SPARSE:%=$(EXAMPLESDIR)/%.o)


ricc: $(RICC_O) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@  -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)

colon: $(COLON_O) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@  -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)

expm: $(EXPM_O) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@  -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)


rand: $(RAND_O) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@  -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)

sparse: $(SPARSE_O) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@  -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)

clean:

	rm -f *.o 

distclean:
	rm -f ricc colon expm rand sparse 


all: ricc colon expm rand sparse

$(EXAMPLESDIR)/%.o : $(EXAMPLESDIR)/%.cxx
	$(CXX) -c $(CXXFLAGS) $< -o $@



