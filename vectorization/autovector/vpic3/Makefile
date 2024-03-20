#Makefile for 3D Electrostatic Vector PIC codes

# Makefile gfortran compiler with MacOS X

#CC = gcc -fopenmp

#CCOPTS = -O3 -Wall -std=c11
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c11
#LEGACY =

# Makefile Intel compiler with Mac OS X

#CC = gcc -fopenmp

#CCOPTS = -O3 -std=c11
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c11
#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

#CC = icc -qopenmp

#CCOPTS = -O3 -std=c11
#CCOPTS = -O3 -mmic -std=c11
#CCOPTS = -O3 -no-vec -std=c11
#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

CC = gcc -fopenmp

CCOPTS = -O3 -ftree-vectorize -Wall -std=c11
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c11
#LEGACY =

# Makefile PGI compiler with Linux
#CC = gcc

#CCOPTS = -O3 -Wall -std=c11
#LEGACY = -Mnomain

# Makefile Cray compiler with Linux

#CC = cc

#CCOPTS = -O 3 -h c11 -h conform
#LEGACY =

#

# Linkage rules

all :  cvpic3

cvpic3 : cvpic3.o cvpush3.o dtimer.o
	$(CC) $(CCOPTS) -o cvpic3 cvpic3.o cvpush3.o \
	dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

cvpush3.o : vpush3.c
	$(CC) $(CCOPTS) -o cvpush3.o -c vpush3.c

cvpic3.o : vpic3.c
	$(CC) $(CCOPTS) -o cvpic3.o -c vpic3.c

clean :
	rm -f *.o *.mod cvpic3

clobber: clean
	rm -f cvpic3 
