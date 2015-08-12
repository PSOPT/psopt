include ../Makefile_linux.inc

STC1 = stc1   $(SNOPT_WRAPPER)

STC1_O = $(STC1:%=$(EXAMPLESDIR)/%.o)


stc1: $(STC1_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

