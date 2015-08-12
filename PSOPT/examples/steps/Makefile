include ../Makefile_linux.inc

STEPS = steps   $(SNOPT_WRAPPER)

STEPS_O = $(STEPS:%=$(EXAMPLESDIR)/%.o)


steps: $(STEPS_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

