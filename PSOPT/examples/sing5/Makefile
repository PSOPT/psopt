include ../Makefile_linux.inc

SING5 = singular5   $(SNOPT_WRAPPER)

SING5_O = $(SING5:%=$(EXAMPLESDIR)/%.o)


sing5: $(SING5_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

