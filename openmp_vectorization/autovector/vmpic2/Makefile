#Makefile for 2D Electrostatic OpenMP/Vector PIC codes

# Makefile gfortran compiler with Linux

# OpenMP
MPFC = gfortran -fopenmp
MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc

FC90 = gfortran
#FC03 = gfortran
CC = gcc

OPTS90 = -O3 -ftree-vectorize
#OPTS90 = -O3 -ftree-vectorizer-verbose=6
#OPTS90 = -O3 -fno-tree-vectorize
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#OPTS03 = -O3 -ftree-vectorize
#OPTS03 = -O3 -ftree-vectorizer-verbose=6
#OPTS03 = -O3 -fno-tree-vectorize
#OPTS03 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS03 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -ftree-vectorize -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

LOPTS = -lgomp

# Linkage rules

all : cvmpic2

cvmpic2 : cvmpic2.o cvmpush2.o complib.o dtimer.o
	$(MPCC) $(CCOPTS) -o cvmpic2 cvmpic2.o cvmpush2.o complib.o \
    dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c


complib.o : omplib.c
	$(MPCC) $(CCOPTS) -o complib.o -c omplib.c


cvmpush2.o : vmpush2.c
	$(MPCC) $(CCOPTS) -o cvmpush2.o -c vmpush2.c


cvmpic2.o : vmpic2.c
	$(CC) $(CCOPTS) -o cvmpic2.o -c vmpic2.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f cvmpic2 
