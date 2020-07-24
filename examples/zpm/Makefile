include ../Makefile_linux.inc

ZPM = zpm   $(SNOPT_WRAPPER)

ZPM_O = $(ZPM:%=$(EXAMPLESDIR)/%.o)


zpm: $(ZPM_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

