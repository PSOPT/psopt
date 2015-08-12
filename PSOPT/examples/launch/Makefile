include ../Makefile_linux.inc

LAUNCH = launch   $(SNOPT_WRAPPER)

LAUNCH_O = $(LAUNCH:%=$(EXAMPLESDIR)/%.o)


launch: $(LAUNCH_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

