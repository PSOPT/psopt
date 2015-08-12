include ../Makefile_linux.inc

SHUTT = shuttle_reentry1    $(SNOPT_WRAPPER)

SHUTT_O = $(SHUTT:%=$(EXAMPLESDIR)/%.o)


shutt: $(SHUTT_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

