include ../Makefile_linux.inc

MPEC = mpec   $(SNOPT_WRAPPER)

MPEC_O = $(MPEC:%=$(EXAMPLESDIR)/%.o)


mpec: $(MPEC_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

