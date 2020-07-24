include ../Makefile_linux.inc

HYPER = hypersensitive   $(SNOPT_WRAPPER)

HYPER_O = $(HYPER:%=$(EXAMPLESDIR)/%.o)


hyper: $(HYPER_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

