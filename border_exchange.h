//
// Created by Magnus Poppe Wang on 27/09/2017.
//

#ifndef TDT4200_PS2_BORDER_EXCHANGE_H
#define TDT4200_PS2_BORDER_EXCHANGE_H

#endif //TDT4200_PS2_BORDER_EXCHANGE_H

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include "RPS_MPI.h"

cell* get_first_row(cell **petri);
cell* get_first_col(cell **petri, int xSize, int ySize);
cell* get_last_row(cell **petri, int ySize);
cell* get_last_col(cell **petri,  int xSize, int ySize);

void exchange_borders(cell** matrix, int xSize, int ySize, int rank, int size);
