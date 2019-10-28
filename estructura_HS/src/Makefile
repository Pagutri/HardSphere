CC = gcc
CFLAGS   = -O3 -Wall -Wstrict-prototypes -ansi -pedantic -lm
EXECUTABLE = ../bin/Sk_HS


$(EXECUTABLE): Sk_HS.o
	$(CC) -o $(EXECUTABLE) Sk_HS.o $(CFLAGS)


Sk_HS.o: Sk_HS.c
	$(CC) -c -o Sk_HS.o Sk_HS.c $(CFLAGS)


clean:
	rm -f *.o $(EXECUTABLE)
