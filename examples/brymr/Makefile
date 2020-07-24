include ../Makefile_linux.inc

BRYMR = bryson_max_range   $(SNOPT_WRAPPER)

BRYMR_O = $(BRYMR:%=$(EXAMPLESDIR)/%.o)


brymr: $(BRYMR_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

