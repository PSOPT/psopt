include ../Makefile_linux.inc

LOWTHR = low_thrust   $(SNOPT_WRAPPER)

LOWTHR_O = $(LOWTHR:%=$(EXAMPLESDIR)/%.o)


lowthr: $(LOWTHR_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

