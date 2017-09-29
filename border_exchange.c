//
// Created by Magnus Poppe Wang on 29/09/2017.
//

#include "border_exchange.h"

int p_north, p_south, p_east, p_west;
MPI_Datatype border_row_t;
MPI_Datatype border_col_t;
MPI_Comm cart_comm;



void send(int* payload, int direction, MPI_Datatype type, int len)
{
    // PRINTING FOR DEBUG:
//    if (direction == p_west) printf("BEFORE SENDING TO WEST:     ");
//    else                     printf("BEFORE SENDING TO EAST:     ");
//    print_list(payload, len);

    // SENDING AND FREEING MEMORY
    MPI_Send(payload, 1, type, direction, 0, cart_comm);
    if (type == border_col_t) free( payload );
}

void recieve(int direction, MPI_Datatype type, int* package, int len)
{
    // RECIEVING:
    MPI_Recv(package, 1, type, direction, 0, cart_comm, MPI_STATUS_IGNORE);

//    // PRINT FOR DEBUG:
//    if (direction == p_north) printf("AFTER RECIEVING FROM NORTH: ");
//    else                      printf("AFTER RECIEVING FROM SOUTH: ");
//    print_list(package, len);
}

void send_row(int row, int direction, int** matrix, int len)
{
    int* payload = matrix[row];
    send(payload, direction, border_row_t, len);
}

void send_col(int col, int direction, int** matrix, int len)
{
    // GETTING THE CORRECT VALUES FOR THE COLUMN
    int* payload = malloc(sizeof(int)*len);
    for (int y = 0; y < len; y++) payload[y] = matrix[y][col];
    send(payload, direction, border_col_t, len);
}

void exchange_borders(
        int **matrix, int xSize, int ySize, int rank, int size,  int n, int s, int e, int w,
        MPI_Datatype row, MPI_Datatype col, MPI_Comm communicator
)
{
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
        send_row(1, p_south, matrix, xSize);

        // RECIEVE PACKAGE FROM SOUTH ONLY:
        int *recieve_south = malloc(xSize * sizeof(int));
        recieve(p_south, border_row_t, recieve_south, xSize);
        stitch_bottom_row(recieve_south, ySize, xSize, matrix);
    }
    else if ((squared * (squared - 1)) + 1 <= rank + 1 && rank + 1 <= size)
    {
        // SEND PACKAGE NORTH ONLY:
        send_row(xSize-2, p_north, matrix, xSize);

        // RECIEVE PACKAGE FROM NORTH ONLY:
        int *recieve_north = malloc(xSize * sizeof(int));
        recieve(p_north, border_row_t, recieve_north, xSize);
        stitch_top_row(recieve_north, xSize, matrix);
    }
    else
    {
        // SEND PACKAGE BOTH NORTH AND SOUTH ONLY:
        send_row(1,       p_south, matrix, xSize);
        send_row(xSize-2, p_north, matrix, xSize);

        // RECIEVE PACKAGE BOTH NORTH AND SOUTH ONLY:
        int *recieve_south = malloc(xSize * sizeof(int));
        recieve(p_south, border_row_t, recieve_south, xSize);
        stitch_bottom_row(recieve_south, ySize, xSize, matrix);

        int *recieve_north = malloc(xSize * sizeof(int));
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
            int* package = malloc(sizeof(int)*ySize);
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
            int* package = malloc(sizeof(int)*ySize);
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
        int* west_package = malloc(sizeof(int)*ySize);
        int* east_package = malloc(sizeof(int)*ySize);
        recieve(p_west, border_col_t, west_package, ySize);
        recieve(p_east, border_col_t, east_package, ySize);

        stitch_left_column(west_package, ySize, xSize, matrix);
        stitch_right_column(east_package, ySize, xSize, matrix);
    }
}

void stitch_bottom_row(int *row, int ylen, int xlen, int** matrix) {
    for (int i = 0; i < xlen; i++) {
        matrix[ylen - 1][i] = row[i];
    }
}

void stitch_top_row(int *row, int len, int** matrix) {
    for (int i = 0; i < len; i++) {
        matrix[0][i] = row[i];
    }
}

void stitch_left_column(int *column, int ylen, int xlen, int** matrix) {
    for (int i = 0; i < ylen; i++) {
        matrix[i][0] = column[i];
    }
}

void stitch_right_column(int *column, int ylen, int xlen, int** matrix) {
    for (int i = 0; i < ylen; i++) {
        matrix[i][xlen-1] = column[i];
    }
}



void print_matrix(int **matrix, int Xlen, int Ylen) {
    for (int y = 0; y < Ylen; y++) {
        for (int x = 0; x < Xlen; x++) {
            if (matrix[y][x] < 10) printf("0%d ", matrix[y][x]);
            else printf("%d ", matrix[y][x]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_list(int *list, int len) {
    for (int i = 0; i < len; i++) {
        if (list[i] < 10) printf("0%d ", list[i]);
        else              printf("%d ",  list[i]);
    }
    printf("\n");
}