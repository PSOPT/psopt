include ../Makefile_linux.inc

TWOPHSC = twophase_schwartz   $(SNOPT_WRAPPER)

TWOPHSC_O = $(TWOPHSC:%=$(EXAMPLESDIR)/%.o)


twophsc: $(TWOPHSC_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

