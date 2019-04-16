include ../Makefile_linux.inc

USER = user   $(SNOPT_WRAPPER)

USER_O = $(USER:%=$(EXAMPLESDIR)/%.o)


user: $(USER_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

