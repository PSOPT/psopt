include ../Makefile_linux.inc

LAMBERT = lambert   $(SNOPT_WRAPPER)

LAMBERT_O = $(LAMBERT:%=$(EXAMPLESDIR)/%.o)


lambert: $(LAMBERT_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

