# CFLAGS+= -std=c99 -g -O3

util_objects = help.o
program = border_exchange.o

SHELL := /bin/bash

all : $(util_objects) $(program)
	mpicc $(util_objects) $(program) -o border_exchange

.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch

clear :
	clear

.PHONY : run
run : clean all clear
	mpirun -n 4 ./border_exchange