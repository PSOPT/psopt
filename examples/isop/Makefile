include ../Makefile_linux.inc

ISOP = isoperimetric   $(SNOPT_WRAPPER)

ISOP_O = $(ISOP:%=$(EXAMPLESDIR)/%.o)


isop: $(ISOP_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

