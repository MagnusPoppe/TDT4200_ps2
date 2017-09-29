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
        int **matrix, int xSize, int ySize, int rank, int size,  int n, int s, int e, int w,
        MPI_Datatype row, MPI_Datatype col, MPI_Comm communicator
);

// FUNCTIONS FOR COMMUNICATION
void send(int* payload, int direction, MPI_Datatype type, int len);
void recieve(int direction, MPI_Datatype type, int* package, int len);
void send_row(int row, int direction, int** matrix, int len);
void send_col(int col, int direction, int** matrix, int len);

// FUNCTIONS FOR STITCHING THE IMAGE TOGETHER:
void stitch_bottom_row(int *row, int ylen, int xlen, int** matrix);
void stitch_top_row(int *row, int len, int** matrix);
void stitch_left_column(int *column, int len, int** matrix);
void stitch_right_column(int *column, int ylen, int xlen, int** matrix);

// UTILS:
void print_list(int *list, int len);
