SHELL         = /bin/sh
USERHOME      = .
F2CINC = $(F2CINCLUDE)
F2CLIB = $(F2CLIBRARY)
DMATRIXDIR=./dmatrix
CXSPARSE=./SuiteSparse/CXSparse
LUSOL=./lusol/csrc

prefix = $(USERHOME)/Ipopt
# Directory with header files
IPOPTINCDIR = ${prefix}/include/coin
# Directory with libipopt.a
IPOPTLIBDIR = ${exec_prefix}/lib
exec_prefix = ${prefix}

PSOPTDIR    = ./PSOPT

PSOPTSRCDIR   = $(PSOPTDIR)/src
EXAMPLESDIR   = $(PSOPTDIR)/examples
CXSPARSE_LIBS = $(CXSPARSE)/Lib/libcxsparse.a
LUSOL_LIBS    = $(LUSOL)/liblusol.a
SPARSE_LIBS   = $(CXSPARSE_LIBS) $(LUSOL_LIBS) -ldl
ADOLC_LIBS    = -ladolc
PSOPT_LIBS    = $(PSOPTDIR)/lib/libpsopt.a
DMATRIX_LIBS  = $(DMATRIXDIR)/lib/libdmatrix.a


$(CXSPARSE_LIBS):
# 	(cp UFconfig.h $(CXSPARSE)/Include;cd $(CXSPARSE)/Lib; $(MAKE))
	(cd $(CXSPARSE);cd ..;$(MAKE))

$(LUSOL_LIBS):
	(cp Makefile.lusol $(LUSOL)/Makefile; cd $(LUSOL); $(MAKE))

$(DMATRIX_LIBS): $(DMATRIXDIR)/src/dmatrixv.cxx $(DMATRIXDIR)/include/dmatrixv.h
	(cd $(DMATRIXDIR)/lib; $(MAKE))


dmatrix_examples:
	(cd $(DMATRIXDIR)/examples; $(MAKE) all)

$(PSOPT_LIBS):
	(cd $(PSOPTDIR)/lib; $(MAKE))

bioreactor:
	(cd $(EXAMPLESDIR)/$@; make $@)

brac1:
	(cd $(EXAMPLESDIR)/$@; make $@)


shutt:
	(cd $(EXAMPLESDIR)/$@; make $@)


manutec:
	(cd $(EXAMPLESDIR)/$@; make $@)

missile:
	(cd $(EXAMPLESDIR)/$@; make $@)

moon:
	(cd $(EXAMPLESDIR)/$@; make $@)

stc1:
	(cd $(EXAMPLESDIR)/$@; make $@)


brymr:
	(cd $(EXAMPLESDIR)/$@; make $@)

twoburn:
	(cd $(EXAMPLESDIR)/$@; make $@)

twolink:
	(cd $(EXAMPLESDIR)/$@; make $@)

twophsc:
	(cd $(EXAMPLESDIR)/$@; make $@)

twophro:
	(cd $(EXAMPLESDIR)/$@; make $@)

hyper:
	(cd $(EXAMPLESDIR)/$@; make $@)


launch:
	(cd $(EXAMPLESDIR)/$@; make $@)

lambert:
	(cd $(EXAMPLESDIR)/$@; make $@)

bryden:
	(cd $(EXAMPLESDIR)/$@; make $@)

delay1:
	(cd $(EXAMPLESDIR)/$@; make $@)

goddard:
	(cd $(EXAMPLESDIR)/$@; make $@)

steps:
	(cd $(EXAMPLESDIR)/$@; make $@)

sing5:
	(cd $(EXAMPLESDIR)/$@; make $@)

climb:
	(cd $(EXAMPLESDIR)/$@; make $@)

cracking:
	(cd $(EXAMPLESDIR)/$@; make $@)

isop:
	(cd $(EXAMPLESDIR)/$@; make $@)

catmix:
	(cd $(EXAMPLESDIR)/$@; make $@)

chain:
	(cd $(EXAMPLESDIR)/$@; make $@)

obstacle:
	(cd $(EXAMPLESDIR)/$@; make $@)

crane:
	(cd $(EXAMPLESDIR)/$@; make $@)

ipc:
	(cd $(EXAMPLESDIR)/$@; make $@)

alpine:
	(cd $(EXAMPLESDIR)/$@; make $@)

lts:
	(cd $(EXAMPLESDIR)/$@; make $@)

user:
	(cd $(EXAMPLESDIR)/$@; make $@)

coulomb:
	(cd $(EXAMPLESDIR)/$@; make $@)

lowthr:
	(cd $(EXAMPLESDIR)/$@; make $@)

heat:
	(cd $(EXAMPLESDIR)/$@; make $@)

zpm:
	(cd $(EXAMPLESDIR)/$@; make $@)

glider:
	(cd $(EXAMPLESDIR)/$@; make $@)

notorious:
	(cd $(EXAMPLESDIR)/$@; make $@)

reorientation:
	(cd $(EXAMPLESDIR)/$@; make $@)

mpec:
	(cd $(EXAMPLESDIR)/$@; make $@)

dae_i3:
	(cd $(EXAMPLESDIR)/$@; make $@)

breakwell:
	(cd $(EXAMPLESDIR)/$@; make $@)

rayleigh:
	(cd $(EXAMPLESDIR)/$@; make $@)

test: launch
	(cd $(EXAMPLESDIR)/launch; ./launch)


all: $(CXSPARSE_LIBS) $(DMATRIX_LIBS) $(LUSOL_LIBS) $(PSOPT_LIBS) dmatrix_examples bioreactor brac1 shutt manutec missile moon stc1 sing5 steps brymr twoburn twolink twophsc twophro hyper launch lambert bryden delay1 goddard sing5 climb cracking isop catmix chain obstacle crane ipc alpine lts user  coulomb lowthr heat zpm glider notorious reorientation mpec dae_i3 breakwell rayleigh test


clean:

	(cd $(DMATRIXDIR)/lib; $(MAKE) clean)
	(cd $(CXSPARSE)/Lib; $(MAKE) clean)
	(cd $(LUSOL); $(MAKE) clean)
	(cd $(DMATRIXDIR)/examples; $(MAKE) clean)
	(cd $(PSOPTDIR)/lib; $(MAKE) clean)


distclean:

	(cd $(DMATRIXDIR)/lib; $(MAKE) distclean)
	(cd $(CXSPARSE)/Lib; $(MAKE) distclean)
	(cd $(LUSOL); $(MAKE) clean)
	(cd $(DMATRIXDIR)/examples; $(MAKE) distclean)
	(cd $(PSOPTDIR)/lib; $(MAKE) distclean)


