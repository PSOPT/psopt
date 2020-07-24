include ../Makefile_linux.inc

REORIENTATION = reorientation   $(SNOPT_WRAPPER)

REORIENTATION_O = $(REORIENTATION:%=$(EXAMPLESDIR)/%.o)


reorientation: $(REORIENTATION_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

