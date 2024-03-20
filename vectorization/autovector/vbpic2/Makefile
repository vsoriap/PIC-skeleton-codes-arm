#Makefile for 2-1/2D Electromagnetic Vector PIC codes

# Makefile gcc compiler with MacOS X

#CC = gcc

#CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

# Makefile Intel compiler with Mac OS X

#CC = gcc

#CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

# Makefile Intel compiler with Linux

#CC = icc

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -std=c99

# Makefile gfortran compiler with Linux

CC = gcc -fopenmp

CCOPTS = -O3 -ftree-vectorize -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

# Makefile PGI compiler with Linux

#CC = gcc

#CCOPTS = -O3 -Wall -std=c99

# Makefile Cray compiler with Linux

#CC = cc

#CCOPTS = -O 3 -h c99 -h conform

#

# Linkage rules

all : cvbpic2

cvbpic2 : cvbpic2.o cvbpush2.o dtimer.o
	$(CC) $(CCOPTS) -o cvbpic2 cvbpic2.o cvbpush2.o \
	dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

cvbpush2.o : vbpush2.c
	$(CC) $(CCOPTS) -o cvbpush2.o -c vbpush2.c

cvbpic2.o : vbpic2.c
	$(CC) $(CCOPTS) -o cvbpic2.o -c vbpic2.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fvbpic2 cvbpic2 f03vbpic2
