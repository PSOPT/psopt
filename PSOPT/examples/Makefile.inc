CXX = cl

CC = cl

LD  = link

CFLAGS = -DWIN32 -DLAPACK -DSPARSE_MATRIX -DUSE_IPOPT  -D_CRT_ERRNO_DEFINED

BASE = ..\..\..

DMATRIX  = $(BASE)\dmatrix

LUSOL    = $(BASE)\lusol\csrc

CXSPARSE = $(BASE)\CXSparse


ADOLC    = C:\ADOL-C-2.1.12\ADOL-C\

ADOLCLIB = C:\ADOL-C-2.1.12\windows\sparse\adolc.lib

IPOPT = C:\Ipopt-3.9.3

include ..\..\..\MicrosoftVisualStudio.inc

PSOPT = $(BASE)\PSOPT

INCLUDE = -I$(PSOPT)\src -I$(DMATRIX)\include  -I$(BASE)\CXSparse\Include -I$(LUSOL) -I$(IPOPT)\include\coin -I$(IPOPT)\Ipopt\src\Interfaces -I$(IPOPT)\Ipopt\src\Common -I$(IPOPT)\BuildTools\headers -I$(IPOPT)\Ipopt\src\LinAlg -I$(ADOLC) -I$(ADOLC)\adolc -I$(ADOLC)\adolc\drivers -I$(ADOLC)\adolc\sparse -I$(VCINCLUDE) -I$(VCSDK)



LIBS = $(DMATRIX)\lib\libdmatrix.lib $(PSOPT)\lib\libpsopt.lib $(DMATRIX)\lib\clapack.lib $(DMATRIX)\lib\blas.lib $(DMATRIX)\lib\libf2c.lib $(IPOPT)\Ipopt\MSVisualStudio\v8\libIpopt\Release\libIpopt.lib $(IPOPT)\Ipopt\MSVisualStudio\v8\libf2c\vcf2c.lib  $(IPOPT)\Ipopt\MSVisualStudio\v8\LibCoinLapack\Release\libCoinLapack.lib  $(IPOPT)\Ipopt\MSVisualStudio\v8\libCoinBlas\Release\libCoinBlas.lib $(IPOPT)\Ipopt\MSVisualStudio\v8\libCoinHSL\Release\libCoinHSL.lib $(LUSOL)\liblusol.lib $(CXSPARSE)\Lib\libcxsparse.lib $(ADOLCLIB)


SNFW_SRC =

SNFW_OBJ =



.cxx.obj:
	$(CXX) -c $(CFLAGS) $(INCLUDE) %s /Fo$*.obj

.c.obj:
	$(CC) -c $(CFLAGS) $(INCLUDE) %s /Fo$*.obj


clean:
	del  *.obj






