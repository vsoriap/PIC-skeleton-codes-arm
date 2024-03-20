#Makefile for 3D Electromagnetic Vector PIC codes

# Makefile gcc compiler with MacOS X

#CC = gcc

#CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99
LEGACY =

# Makefile Intel compiler with Mac OS X

#CC = gcc

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99
#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

#CC = icc

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -mmic -std=c99
#CCOPTS = -O3 -no-vec -std=c99
#LEGACY = -nofor_main

# Makefile gcc compiler with Linux

CC = gcc -fopenmp

CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99
#LEGACY =

# Makefile PGI compiler with Linux

#CC = gcc -fopenmp

#CCOPTS = -O3 -Wall -std=c99
#LEGACY = -Mnomain

# Makefile Cray compiler with Linux

#CC = cc

#CCOPTS = -O 3 -h c99 -h conform
#LEGACY =

#

# Linkage rules

#all : fvbpic3 cvbpic3 f03vbpic3
all : cvbpic3


cvbpic3 : cvbpic3.o cvbpush3.o dtimer.o
	$(CC) $(CCOPTS) -o cvbpic3 cvbpic3.o cvbpush3.o \
	dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

cvbpush3.o : vbpush3.c
	$(CC) $(CCOPTS) -o cvbpush3.o -c vbpush3.c

cvbpic3.o : vbpic3.c
	$(CC) $(CCOPTS) -o cvbpic3.o -c vbpic3.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f cvbpic3