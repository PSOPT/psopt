include ../Makefile_linux.inc

IPC = interior_point   $(SNOPT_WRAPPER)

IPC_O = $(IPC:%=$(EXAMPLESDIR)/%.o)


ipc: $(IPC_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)
	rm -f *.o

