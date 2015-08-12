include ../Makefile_linux.inc

DELAY1 = delay1   $(SNOPT_WRAPPER)

DELAY1_O = $(DELAY1:%=$(EXAMPLESDIR)/%.o)


delay1: $(DELAY1_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

