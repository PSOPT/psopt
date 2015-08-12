include ../Makefile_linux.inc

TWOPHRO = twophase_robot   $(SNOPT_WRAPPER)

TWOPHRO_O = $(TWOPHRO:%=$(EXAMPLESDIR)/%.o)


twophro: $(TWOPHRO_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

