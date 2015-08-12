include ../Makefile_linux.inc

MISSILE = missile   $(SNOPT_WRAPPER)

MISSILE_O = $(MISSILE:%=$(EXAMPLESDIR)/%.o)


missile: $(MISSILE_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

