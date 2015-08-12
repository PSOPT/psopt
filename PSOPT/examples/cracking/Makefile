include ../Makefile_linux.inc

CRACKING = cracking   $(SNOPT_WRAPPER)

CRACKING_O = $(CRACKING:%=$(EXAMPLESDIR)/%.o)


cracking: $(CRACKING_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

