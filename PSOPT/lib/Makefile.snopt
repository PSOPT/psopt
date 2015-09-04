
#Usage:  make
#add ../lib to LD_LIBRARY_PATH

SHELL         = /bin/sh

USERHOME      = /home/$(shell whoami)

F2CINC = $(F2CINCLUDE)
F2CLIB = $(F2CLIBRARY)

DMATRIXDIR=../../dmatrix
CXSPARSE=../../SuiteSparse/CXSparse
LUSOL=../../lusol/csrc

SNOPTDIR=$(USERHOME)/snopt7
SNOPTINTERFACEDIR=../../snopt-interface

prefix = $(USERHOME)/Ipopt-3.12.3
# Directory with header files
IPOPTINCDIR = -I${prefix}/Ipopt/src/Interfaces -I${prefix}/Ipopt/src/Common -I${prefix}/Ipopt/src/LinAlg
# Directory with libipopt.a
IPOPTLIBDIR = ${exec_prefix}/lib
exec_prefix = ${prefix}


PSOPTSRCDIR = ../src


CXX           = /usr/bin/g++
CC            = /usr/bin/gcc
CXXFLAGS      = -O0 -g -I$(SNOPTINTERFACEDIR)/include -I$(USERHOME)/adolc_base/include  -I$(DMATRIXDIR)/include -I$(PSOPTSRCDIR) -DLAPACK -DUNIX -DSPARSE_MATRIX -DUSE_SNOPT -DUSE_IPOPT -I$(CXSPARSE)/Include -I$(CXSPARSE)/../SuiteSparse_config -I$(LUSOL) $(IPOPTINCDIR) -fomit-frame-pointer -pipe -DNDEBUG -pedantic-errors -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -fPIC -DHAVE_MALLOC

CFLAGS        = -O0 -fPIC


PSOPTLIB = libpsopt.a

$(PSOPTLIB):  $(PSOPTLIB)($(PSOPTSRCDIR)/psopt.o $(PSOPTSRCDIR)/plot.o $(PSOPTSRCDIR)/util.o $(PSOPTSRCDIR)/pseudospectral.o $(PSOPTSRCDIR)/propagate.o $(PSOPTSRCDIR)/print.o $(PSOPTSRCDIR)/validate.o $(PSOPTSRCDIR)/scaling.o $(PSOPTSRCDIR)/interpolation.o $(PSOPTSRCDIR)/NLP_objective.o $(PSOPTSRCDIR)/NLP_constraints.o $(PSOPTSRCDIR)/mesh.o $(PSOPTSRCDIR)/evaluate.o $(PSOPTSRCDIR)/workspace.o $(PSOPTSRCDIR)/get_numbers.o $(PSOPTSRCDIR)/get_variables.o $(PSOPTSRCDIR)/setup.o $(PSOPTSRCDIR)/solution.o $(PSOPTSRCDIR)/NLP_guess.o $(PSOPTSRCDIR)/NLP_bounds.o $(PSOPTSRCDIR)/NLP_interface.o  $(PSOPTSRCDIR)/IPOPT_interface.o $(PSOPTSRCDIR)/derivatives.o $(PSOPTSRCDIR)/trajectories.o $(PSOPTSRCDIR)/SNOPT_interface.o $(PSOPTSRCDIR)/user_functions.o $(PSOPTSRCDIR)/integrate.o $(PSOPTSRCDIR)/phases.o $(PSOPTSRCDIR)/parameter_estimation.o)


clean:
	rm -f *.o

distclean:
	rm -f *.a


$(PSOPTSRCDIR)/%.o : $(PSOPTSRCDIR)/%.cxx
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(DMATRIXDIR)/src/%.o : $(DMATRIXDIR)/src/%.cxx
	$(CXX) -c $(CXXFLAGS) $< -o $@


