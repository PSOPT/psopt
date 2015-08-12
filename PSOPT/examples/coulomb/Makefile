include ../Makefile_linux.inc

COULOMB = coulomb   $(SNOPT_WRAPPER)

COULOMB_O = $(COULOMB:%=$(EXAMPLESDIR)/%.o)


coulomb: $(COULOMB_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o	

