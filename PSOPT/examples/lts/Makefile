include ../Makefile_linux.inc

LTS = lts   $(SNOPT_WRAPPER)

LTS_O = $(LTS:%=$(EXAMPLESDIR)/%.o)


lts: $(LTS_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

