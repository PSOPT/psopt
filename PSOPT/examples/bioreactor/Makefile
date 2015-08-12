include ../Makefile_linux.inc

BIOREACTOR = bioreactor   $(SNOPT_WRAPPER)

BIOREACTOR_O = $(BIOREACTOR:%=$(EXAMPLESDIR)/%.o)


bioreactor: $(BIOREACTOR_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

