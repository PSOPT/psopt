# Modify the "-O" optimization option for best performance (-O3 on Linux):
CC = gcc
CFLAGS = -DYZHANG -DNDEBUG -Wall -fexceptions -pthread -ansi -fPIC -O3
I = -I.

AR = ar cr
RANLIB = ranlib

all: liblusol.a 

LUSOL_SOURCE  =  lusol.c  mmio.c commonlib.c lusolio.c hbio.c myblas.c 
LUSOL_OBJ     =  lusol.o  mmio.o commonlib.o lusolio.o hbio.o myblas.o 
LUSOL_LIB_OBJ =  lusol.o  mmio.o commonlib.o lusolio.o hbio.o myblas.o 

LUSOL_LIB    = liblusol.a


.c.o: %.c
	$(CC) $(CFLAGS) $(I) -c $< -o $@

$(LUSOL_LIB): $(LUSOL_LIB_OBJ)
	$(AR) $(LUSOL_LIB) $(LUSOL_LIB_OBJ)
	$(RANLIB) $(LUSOL_LIB)
	cp liblusol.a ..


clean:
	rm -f *.o

purge: distclean

distclean: 
	rm -f *.a
