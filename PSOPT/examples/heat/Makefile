include ../Makefile_linux.inc

HEAT = heat   $(SNOPT_WRAPPER)

HEAT_O = $(HEAT:%=$(EXAMPLESDIR)/%.o)


heat: $(HEAT_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

