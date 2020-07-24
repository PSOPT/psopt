include ../Makefile_linux.inc

TWOLINK = twolinkarm   $(SNOPT_WRAPPER)

TWOLINK_O = $(TWOLINK:%=$(EXAMPLESDIR)/%.o)


twolink: $(TWOLINK_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

