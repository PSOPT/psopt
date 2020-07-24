include ../Makefile_linux.inc

ALPINE = alpine   $(SNOPT_WRAPPER)

ALPINE_O = $(ALPINE:%=$(EXAMPLESDIR)/%.o)


alpine: $(ALPINE_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

