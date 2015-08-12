include ../Makefile_linux.inc

CLIMB = climb   $(SNOPT_WRAPPER)

CLIMB_O = $(CLIMB:%=$(EXAMPLESDIR)/%.o)


climb: $(CLIMB_O)  $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

