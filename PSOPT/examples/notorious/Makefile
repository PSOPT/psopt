include ../Makefile_linux.inc

NOTORIOUS = notorious   $(SNOPT_WRAPPER)

NOTORIOUS_O = $(NOTORIOUS:%=$(EXAMPLESDIR)/%.o)


notorious: $(NOTORIOUS_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

