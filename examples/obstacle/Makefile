include ../Makefile_linux.inc

OBSTACLE = obstacle   $(SNOPT_WRAPPER)

OBSTACLE_O = $(OBSTACLE:%=$(EXAMPLESDIR)/%.o)


obstacle: $(OBSTACLE_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

