#ifndef RPS_MPI_H
#define RPS_MPI_H


#include "CA.h"
#include "bitmap.h"
#include "RPS.h"
#include "border_exchange.h"

#define WHITE   0
#define ROCK    1
#define PAPER   2
#define SCISSOR 3

#define IMG_X 512
#define IMG_Y 512

// How many iterations?
#define ITERATIONS 10000

#endif
void create_cartesian_communicator();
void initialize();
void iterate_CA(cell** current_image, cell** next_image);
void gather_petri();
void create_types();
void convert_1D_to_2D(cell* whole_petri, int len);
void write_images();
void free_multi_cell_array(cell **array, int width);
