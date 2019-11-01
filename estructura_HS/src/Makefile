CC = gcc
CFLAGS   = -O3 -Wall -Wstrict-prototypes -ansi -pedantic -lm
COMMONSRC = libraries.h constants.h
EXECUTABLE = ../bin/forces_HS


$(EXECUTABLE): forces_HS.o funciones.o $(COMMONSRC)
	$(CC) -o $(EXECUTABLE) forces_HS.o funciones.o $(CFLAGS)

funciones.o: funciones.c $(COMMONSRC)
	$(CC) -c -o funciones.o funciones.c $(CFLAGS)

forces_HS.o: forces_HS.c $(COMMONSRC)
	$(CC) -c -o forces_HS.o forces_HS.c $(CFLAGS)

Sk_HS.o: Sk_HS.c $(COMMONSRC)
	$(CC) -c -o Sk_HS.o Sk_HS.c $(CFLAGS)


clean:
	rm -f *.o $(EXECUTABLE)
