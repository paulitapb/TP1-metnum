CC = g++
CFLAGS =  -pedantic -O3 -lm -g
#-Wall -Wextra
TARGET = main
.PHONY: all
all: tp1

valgrind: 
	valgrind --leak-check=full  --undef-value-errors=yes ./main

all: tp1

tp1: main_cpp.o grafo_cpp.o matriz_rala_cpp.o matriz_ralaCSR_cpp.o tests_matriz_rala_cpp.o utilidades.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main_cpp.o: main.cpp 
	$(CC) $(CFLAGS) -c $< -o $@
grafo_cpp.o: grafo.cpp grafo.h
	$(CC) $(CFLAGS) -c $< -o $@
matriz_rala_cpp.o: matriz_rala.cpp matriz_rala.h utilidades.h
	$(CC) $(CFLAGS) -c $< -o $@
matriz_ralaCSR_cpp.o: matriz_ralaCSR.cpp matriz_ralaCSR.h 
	$(CC) $(CFLAGS) -c $< -o $@
tests_matriz_rala_cpp.o: tests_matriz_rala.cpp tests_matriz_rala.h
	$(CC) $(CFLAGS) -c $< -o $@
utilidades.o: utilidades.cpp utilidades.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm -f *.o
	rm -f tp1