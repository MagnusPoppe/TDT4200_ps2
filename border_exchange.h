//
// Created by Magnus Poppe Wang on 29/09/2017.
//

#ifndef TDT4200_PS2_BORDER_EXCHANGE_H
#define TDT4200_PS2_BORDER_EXCHANGE_H

#endif //TDT4200_PS2_BORDER_EXCHANGE_H

#include "RPS_MPI.h"
#include <mpi.h>

// FUNCTION FOR THE MAIN BORDER EXCHANGE
void exchange_borders(
        cell **matrix, int xSize, int ySize, int rank, int size,  int n, int s, int e, int w,
        MPI_Datatype row, MPI_Datatype col, MPI_Comm communicator
);

// FUNCTIONS FOR COMMUNICATION
void send(cell* payload, int direction, MPI_Datatype type, int len);
void recieve(int direction, MPI_Datatype type, cell* package, int len);
void send_row(int row, int direction, cell** matrix, int len);
void send_col(int col, int direction, cell** matrix, int len);

// FUNCTIONS FOR STITCHING THE IMAGE TOGETHER:
void stitch_bottom_row(cell *row, int ylen, int xlen, cell** matrix);
void stitch_top_row(cell *row, int len, cell** matrix);
void stitch_left_column(cell *column, int ylen, int xlen, cell** matrix);
void stitch_right_column(cell *column, int ylen, int xlen, cell** matrix);

// UTILS:
void print_list(cell *list, int len);
void print_matrix(cell **matrix, int Xlen, int Ylen);
