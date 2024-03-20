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

cvbpic2 : cvbpic2.o csselib2.o cssebpush2.o dtimer.o
	$(CC) $(CCOPTS) -o cvbpic2 cvbpic2.o csselib2.o \
	cssebpush2.o dtimer.o -lm


# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

csselib2.o : sselib2.c
	$(CC) $(CCOPTS) -o csselib2.o -c sselib2.c

csseflib2.o : sseflib2.c
	$(CC) $(CCOPTS) -o csseflib2.o -c sseflib2.c

cssebpush2.o : ssebpush2.c
	$(CC) $(CCOPTS) -o cssebpush2.o -c ssebpush2.c

cvbpic2.o : vbpic2.c
	$(CC) $(CCOPTS) -o cvbpic2.o -c vbpic2.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fvbpic2 cvbpic2 f03vbpic2
