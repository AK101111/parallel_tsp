CC		= gcc

CFLAGS		= -g  -Wall

CLIBS		= -lm 

all: main

main: source/main.c ptsm.o source/ptsm.h
	$(CC) $(CFLAGS) -o main source/main.c ptsm.o

ptsm.o: source/ptsm.c source/ptsm.h
	$(CC) $(CFLAGS) -c source/ptsm.c

clean:		
	rm -f *~ *.o a.out core main ptsm.o  