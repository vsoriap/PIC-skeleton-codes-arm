#Makefile for 2D Electrostatic Vector PIC codes

# Makefile gfortran compiler with MacOS X
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
CC = gcc

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

all : cvpic2

cvpic2 : cvpic2.o cneonlib2.o cneonpush2.o dtimer.o
	$(CC) $(CCOPTS) -o cvpic2 cvpic2.o cneonlib2.o cneonpush2.o \
	dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

cneonlib2.o : neonlib2.c
	$(CC) $(CCOPTS) -o cneonlib2.o -c neonlib2.c

cneonpush2.o : neonpush2.c
	$(CC) $(CCOPTS) -o cneonpush2.o -c neonpush2.c

cvpic2.o : vpic2.c
	$(CC) $(CCOPTS) -o cvpic2.o -c vpic2.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f cvpic2
