include ../Makefile_linux.inc

GODDARD = goddard   $(SNOPT_WRAPPER)

GODDARD_O = $(GODDARD:%=$(EXAMPLESDIR)/%.o)


goddard: $(GODDARD_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

