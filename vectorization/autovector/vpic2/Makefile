#Makefile for 2D Electrostatic Vector PIC codes

# Makefile gcc compiler with MacOS X

#CC = gcc -fopenmp

#CCOPTS = -O3 -Wall -std=c99 
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

# Makefile Intel compiler with Mac OS X

#CC = gcc

#CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

# Makefile Intel compiler with Linux

#CC = icc -qopenmp

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -std=c99

# Makefile gcc compiler with Linux

CC = gcc -fopenmp

CCOPTS = -O3 -ftree-vectorize -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

# Makefile gcc compiler with Linux

#CC = gcc

#CCOPTS = -O3 -Wall -std=c99

# Makefile Cray compiler with Linux

#CC = cc

#OPTS03 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform

#

# Linkage rules

all : cvpic2

cvpic2 : cvpic2.o cvpush2.o dtimer.o
	$(CC) $(CCOPTS) -o cvpic2 cvpic2.o cvpush2.o dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

cvpush2.o : vpush2.c
	$(CC) $(CCOPTS) -o cvpush2.o -c vpush2.c

cvpic2.o : vpic2.c
	$(CC) $(CCOPTS) -o cvpic2.o -c vpic2.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f cvpic2 
