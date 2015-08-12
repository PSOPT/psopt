include ../Makefile_linux.inc

BREAKWELL = breakwell   $(SNOPT_WRAPPER)

BREAKWELL_O = $(BREAKWELL:%=$(EXAMPLESDIR)/%.o)


breakwell: $(BREAKWELL_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

