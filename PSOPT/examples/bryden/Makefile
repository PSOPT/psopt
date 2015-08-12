include ../Makefile_linux.inc

BRYDEN = bryson_denham   $(SNOPT_WRAPPER)

BRYDEN_O = $(BRYDEN:%=$(EXAMPLESDIR)/%.o)


bryden: $(BRYDEN_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

