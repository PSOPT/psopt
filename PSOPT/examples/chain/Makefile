include ../Makefile_linux.inc

CHAIN = chain   $(SNOPT_WRAPPER)

CHAIN_O = $(CHAIN:%=$(EXAMPLESDIR)/%.o)


chain: $(CHAIN_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o
