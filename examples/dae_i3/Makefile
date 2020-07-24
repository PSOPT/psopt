include ../Makefile_linux.inc

DAE_I3 = dae_i3   $(SNOPT_WRAPPER)

DAE_I3_O = $(DAE_I3:%=$(EXAMPLESDIR)/%.o)


dae_i3: $(DAE_I3_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

