include ../Makefile_linux.inc

TWOBURN = twoburn   $(SNOPT_WRAPPER)

TWOBURN_O = $(TWOBURN:%=$(EXAMPLESDIR)/%.o)


twoburn: $(TWOBURN_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

