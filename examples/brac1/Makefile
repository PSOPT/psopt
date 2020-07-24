include ../Makefile_linux.inc

BRAC1 = brac1   $(SNOPT_WRAPPER)

BRAC1_O = $(BRAC1:%=$(EXAMPLESDIR)/%.o)


brac1: $(BRAC1_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

