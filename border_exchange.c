//
// Created by Magnus Poppe Wang on 29/09/2017.
//

#include "border_exchange.h"

int p_north, p_south, p_east, p_west;
MPI_Datatype border_row_t;
MPI_Datatype border_col_t;
MPI_Comm cart_comm;
int rank, size;


void send(cell* payload, int direction, MPI_Datatype type, int len)
{
//     PRINTING FOR DEBUG:
//    if (direction == p_north)      printf("RANK %d, BEFORE SENDING TO NORTH:     ", rank);
//    else if (direction == p_south) printf("RANK %d, BEFORE SENDING TO SOUTH:     ", rank);
//    else if (direction == p_west)  printf("RANK %d, BEFORE SENDING TO WEST:      ", rank);
//    else                           printf("RANK %d, BEFORE SENDING TO EAST:      ", rank);

//    print_list(payload, len);

    // SENDING AND FREEING MEMORY
    MPI_Send(payload, 1, type, direction, 0, cart_comm);
    if (type == border_col_t) free( payload );
}

void recieve(int direction, MPI_Datatype type, cell* package, int len)
{
    // RECIEVING:
    MPI_Recv(package, 1, type, direction, 0, cart_comm, MPI_STATUS_IGNORE);

    // PRINT FOR DEBUG:
//    if (direction == p_north)      printf("RANK %d, AFTER RECIEVING TO NORTH:    ", rank);
//    else if (direction == p_south) printf("RANK %d, AFTER RECIEVING TO SOUTH:    ", rank);
//    else if (direction == p_west)  printf("RANK %d, AFTER RECIEVING TO WEST:     ", rank);
//    else                           printf("RANK %d, AFTER RECIEVING TO EAST:     ", rank);

//    print_list(package, len);
}

void send_row(int row, int direction, cell** matrix, int len)
{
    cell* payload = matrix[row];
    send(payload, direction, border_row_t, len);
}

void send_col(int col, int direction, cell** matrix, int len)
{
    // GETTING THE CORRECT VALUES FOR THE COLUMN
    cell* payload = malloc(sizeof(cell)*len);
    for (int y = 0; y < len; y++)
        payload[y] = matrix[y][col];
    send(payload, direction, border_col_t, len);
}

void exchange_borders(
        cell **matrix, int xSize, int ySize, int rk, int se,       // General info
        int n, int s, int e, int w,                                // Cartesian neighbours
        MPI_Datatype row, MPI_Datatype col, MPI_Comm communicator  // MPI data
)
{
    rank = rk;
    size = se;
    p_north = n;
    p_south = s;
    p_east = e;
    p_west = w;
    border_row_t = row;
    border_col_t = col;
    cart_comm = communicator;

    int squared = sqrt(size);

    if (0 <= rank + 1 && rank + 1 <= squared)
    {
        // SEND PACKAGE SOUTH ONLY:
        send_row(xSize-2, p_south, matrix, xSize);

        // RECIEVE PACKAGE FROM SOUTH ONLY:
        cell *recieve_south = malloc(xSize * sizeof(cell));
        recieve(p_south, border_row_t, recieve_south, xSize);
        stitch_bottom_row(recieve_south, ySize, xSize, matrix);
    }
    else if ((squared * (squared - 1)) + 1 <= rank + 1 && rank + 1 <= size)
    {
        // SEND PACKAGE NORTH ONLY:
        send_row(1, p_north, matrix, xSize);

        // RECIEVE PACKAGE FROM NORTH ONLY:
        cell *recieve_north = malloc(xSize * sizeof(cell));
        recieve(p_north, border_row_t, recieve_north, xSize);
        stitch_top_row(recieve_north, xSize, matrix);
    }
    else
    {
        // SEND PACKAGE BOTH NORTH AND SOUTH ONLY:
        send_row(1,       p_south, matrix, xSize);
        send_row(xSize-2, p_north, matrix, xSize);

        // RECIEVE PACKAGE BOTH NORTH AND SOUTH ONLY:
        cell *recieve_south = malloc(xSize * sizeof(cell));
        recieve(p_south, border_row_t, recieve_south, xSize);
        stitch_bottom_row(recieve_south, ySize, xSize, matrix);

        cell *recieve_north = malloc(xSize * sizeof(cell));
        recieve(p_north, border_row_t, recieve_north, xSize);
        stitch_top_row(recieve_north, xSize, matrix);
    }

    bool not_edge = true;
    for (int i = 0; i < squared; i++)
    {
        if (rank == squared * i)                    // LEFT END
        {
            // SENDING:
            send_col(ySize-2, p_east, matrix, ySize);

            // RECIEVING:
            cell* package = malloc(sizeof(cell)*ySize);
            recieve(p_east, border_col_t, package, ySize);

            // STITCHING PACKAGE TO IMAGE:
            stitch_right_column(package, ySize, xSize, matrix);

            // BREAKING THE LOOP:
            not_edge = false;
            break;
        }
        else if (rank == (squared * (i + 1)) - 1)   // RIGHT END
        {
            // SENDING:
            send_col(1, p_west, matrix, ySize);

            // RECIEVING:
            cell* package = malloc(sizeof(cell)*ySize);
            recieve(p_west, border_col_t, package, ySize);

            // STITCHING PACKAGE TO IMAGE:
            stitch_left_column(package, ySize, xSize, matrix);

            // BREAKING THE LOOP:
            not_edge = false;
            break;
        }
    }
    if (not_edge && squared > 2) // NOT EDGE OF GRID:
    {
        // SENDING:
        send_col(1, p_west, matrix, ySize);
        send_col(ySize-2, p_east, matrix, ySize);

        // RECIEVING:
        cell* west_package = malloc(sizeof(cell)*ySize);
        cell* east_package = malloc(sizeof(cell)*ySize);
        recieve(p_west, border_col_t, west_package, ySize);
        recieve(p_east, border_col_t, east_package, ySize);

        stitch_left_column(west_package, ySize, xSize, matrix);
        stitch_right_column(east_package, ySize, xSize, matrix);
    }
}

void stitch_bottom_row(cell *row, int ylen, int xlen, cell** matrix) {
    for (int i = 0; i < xlen; i++) {
        matrix[ylen - 1] = row;
    }
}

void stitch_top_row(cell *row, int len, cell** matrix) {
    for (int i = 0; i < len; i++) {
        matrix[0] = row;
    }
}

void stitch_left_column(cell *column, int ylen, int xlen, cell** matrix) {
    for (int i = 0; i < ylen; i++) {
        matrix[i][0] = column[i];
    }
}

void stitch_right_column(cell *column, int ylen, int xlen, cell** matrix) {
    for (int i = 0; i < ylen; i++) {
        matrix[i][xlen-1] = column[i];
    }
}



void print_matrix(cell **matrix, int Xlen, int Ylen) {
    for (int y = 0; y < Ylen; y++) {
        for (int x = 0; x < Xlen; x++)
            printf("[%d, %d] ", matrix[y][x].strength, matrix[y][x].color);
        printf("\n");
    }
    printf("\n");
}

void print_list(cell *list, int len) {
    for (int i = 0; i < len; i++)
        printf("[%d, %d] ", list[i].strength, list[i].color);
    printf("\n");
}