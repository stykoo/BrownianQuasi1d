CC=g++
CFLAGS=-W -Wall -ansi -pedantic -std=c++14 -O3
LDFLAGS=-lboost_program_options -lpthread
EXEC=brownianQuasi1d

all: $(EXEC)

$(EXEC): main.o SimulPipe.o simul.o parseArguments.o parameters.o
		$(CC) -o $@ $^ $(LDFLAGS)

simul_test: simul_test.cpp simul.o parameters.o
		$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp parseArguments.h parameters.h
		$(CC) -o $@ -c $< $(CFLAGS)

simul.o: simul.cpp simul.h parameters.h
		$(CC) -o $@ -c $< $(CFLAGS)

SimulPipe.o: SimulPipe.cpp SimulPipe.h simul.h parameters.h
		$(CC) -o $@ -c $< $(CFLAGS)

parseArguments.o: parseArguments.cpp parseArguments.h parameters.h
		$(CC) -o $@ -c $< $(CFLAGS)

parameters.o: parameters.cpp parameters.h
		$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
		rm -f *.o

mrproper: clean
		rm -rf $(EXEC)
