include ../Makefile_linux.inc

RAYLEIGH = rayleigh   $(SNOPT_WRAPPER)

RAYLEIGH_O = $(RAYLEIGH:%=$(EXAMPLESDIR)/%.o)


rayleigh: $(RAYLEIGH_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

