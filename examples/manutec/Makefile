include ../Makefile_linux.inc

MANUTEC = manutec   $(SNOPT_WRAPPER)

MANUTEC_O = $(MANUTEC:%=$(EXAMPLESDIR)/%.o)


manutec: $(MANUTEC_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

