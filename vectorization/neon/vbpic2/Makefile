#Makefile for 2-1/2D Electromagnetic Vector PIC codes

# Makefile gcc compiler with MacOS X

#CC = gcc

#CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

# Makefile gcc compiler with Linux

#CC = gcc

#CCOPTS = -O3 -ftree-vectorize -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

# Linkage rules

all : cvbpic2  

cvbpic2 : cvbpic2.o cneonlib2.o cneonbpush2.o dtimer.o
	$(CC) $(CCOPTS) -o cvbpic2 cvbpic2.o cneonlib2.o \
	cneonbpush2.o dtimer.o -lm


# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

cneonlib2.o : neonlib2.c
	$(CC) $(CCOPTS) -o cneonlib2.o -c neonlib2.c

cneonbpush2.o : neonbpush2.c
	$(CC) $(CCOPTS) -o cneonbpush2.o -c neonbpush2.c

cvbpic2.o : vbpic2.c
	$(CC) $(CCOPTS) -o cvbpic2.o -c vbpic2.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f cvbpic2
