include ../Makefile_linux.inc

CATMIX = catmix   $(SNOPT_WRAPPER)

CATMIX_O = $(CATMIX:%=$(EXAMPLESDIR)/%.o)


catmix: $(CATMIX_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

