include ../Makefile_linux.inc

MOON = moon   $(SNOPT_WRAPPER)

MOON_O = $(MOON:%=$(EXAMPLESDIR)/%.o)


moon: $(MOON_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

