include ../Makefile_linux.inc

GLIDER = glider   $(SNOPT_WRAPPER)

GLIDER_O = $(GLIDER:%=$(EXAMPLESDIR)/%.o)


glider: $(GLIDER_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

