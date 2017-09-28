#ifndef RPS_MPI_H
#define RPS_MPI_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "CA.h"
#include "bitmap.h"
#include "RPS.h"

#define WHITE   0
#define ROCK    1
#define PAPER   2
#define SCISSOR 3

#define IMG_X 4
#define IMG_Y 4

// Each cell is updated based on neighbors of distance 1
#define BORDER_SIZE 1

// How many iterations?
#define ITERATIONS 10000

#endif
