include ../Makefile_linux.inc

CRANE = crane   $(SNOPT_WRAPPER)

CRANE_O = $(CRANE:%=$(EXAMPLESDIR)/%.o)


crane: $(CRANE_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

