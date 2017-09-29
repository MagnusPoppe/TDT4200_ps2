# CFLAGS+= -std=c99 -g -O3

# by uncommenting this line the preprocessor will see #ifdef DEBUG as true
# CFLAGS+= -DDEBUG
util_objects = CA.o bitmap.o
serial_objects = RPS.o
parallel_objects = RPS_MPI.o border_exchange.o

SHELL := /bin/bash

# When running "make" in your terminal, the first target will be chosen.
# The syntax of make is basically:
# name : [stuff I need to make this target]

# In this case the dependencies are easy to figure out, so I will not elaborate further
parallel: $(util_objects) $(parallel_objects)
	mpicc $(util_objects) $(parallel_objects) -o parallelRPS

serial : $(util_objects) $(serial_objects)
	gcc $(util_objects) $(serial_objects) -o serialRPS


RPS_MPI.o: RPS_MPI.c RPS_MPI.h
	mpicc -c RPS_MPI.c RPS_MPI.h

# In this target [stuff I need to make this target] is two other targets, namely clean and all
# This command simply runs the clean target, then the all target, thus recompiling everything.
remake : clean all

# We add .PHONY when a target doesn't actually create any output. In this case we just run a shell
# command, removing all object files, i.e files ending on .o
# the * syntax means [anything].o
.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch && rm -f data/*.bmp && rm -f *.mp4

.PHONY : clear
clear :
	clear

.PHONY : run
run : clean parallel clear
	mpirun -n 4 ./parallelRPS

# Finally, the test target. Builds the 'all' target, then runs the test script on the output
.PHONY : video
video : clean serial
	./serialRPS && ffmpeg -framerate 60 -i data/CA-%000d.bmp -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
