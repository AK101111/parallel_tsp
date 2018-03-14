CC		= clang

CFLAGS		= -g  -Wall

CLIBS		= -lm 

all: main

main: main.c ptsm.o ptsm.h
	$(CC) $(CFLAGS) -o main main.c ptsm.o

ptsm.o: ptsm.c ptsm.h
	$(CC) $(CFLAGS) -c ptsm.c

clean:		
	rm -f *~ *.o a.out core main ptsm.o  