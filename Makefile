CC=g++
CFLAGS=-W -Wall -ansi -pedantic -std=c++14 -O3
LDFLAGS=-lboost_program_options -lpthread
EXEC=brownianQuasi1d
SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
		$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

SimulPipe.o: simul.h

SimulTonks.o: simul.h

%.o: %.cpp %.h parameters.h
		$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
		rm -f *.o

mrproper: clean
		rm -rf $(EXEC)
