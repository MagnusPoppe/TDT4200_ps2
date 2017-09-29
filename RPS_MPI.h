#ifndef RPS_MPI_H
#define RPS_MPI_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "CA.h"
#include "bitmap.h"
#include "RPS.h"
#include "border_exchange.h"

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

void initialize();
//void initialize_petri();
//void iterate_CA();
void gather_petri();
void create_types();
int *get_first_row(int **petri);
void get_first_col(int **petri, int* output, int xSize, int ySize);
int *get_last_row(int **petri, int ySize);
void get_last_col(int **petri, int* output, int xSize, int ySize);
int** create_full_petri(int* whole_petri);
void print_int_list(int *list, int len);