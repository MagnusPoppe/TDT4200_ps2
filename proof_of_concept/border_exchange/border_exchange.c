//
// Created by Magnus Poppe Wang on 27/09/2017.
//



#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "help.h"
#include "../../RPS.h"

#include <mpi.h>

int size;
int rank;
MPI_Comm cart_comm;
int p_north, p_south, p_east, p_west;

int x, y;

int x_dimension, y_dimension;

MPI_Datatype mpi_cell_t;    // Already implemented
void create_cell_type();

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // print_array(array, gridSize);
    // printf("\n\n");

    ////////////////////////////////
    // Create cartesian communicator
    int dims[2];
    dims[0] = x_dimension;
    dims[1] = y_dimension;

    int periods[2]; // we set these to 0 because we are not interested in wrap-around
    periods[0] = 0;
    periods[1] = 0;

    // Local coords:
    int coords[2];
    coords[0] = x;
    coords[1] = y;

    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 2, coords);


    MPI_Cart_shift(cart_comm, 0, 1, &p_north, &p_south);
    MPI_Cart_shift(cart_comm, 1, 1, &p_west, &p_east);

    create_cell_type();


    // A super basic example sending some data:
    int grid = 32 / size;
    int **array = create_2D_array(grid);

//    if (rank == 0) {
//        int recieved;
//        MPI_Recv(&recieved, 1, MPI_INT, p_south, 0, cart_comm, MPI_STATUS_IGNORE);
//        printf("Recieved from the south: %d\n", recieved);
//        recieved++;
//
//        MPI_Send(&recieved, 1, MPI_INT, p_south, 0, cart_comm);
//    }
//    else if (rank == 2) {
//        // Sending 2016 north
//        int sendt = 2016;
//        MPI_Send(&sendt, 1, MPI_INT, p_north, 0, cart_comm);
//
//        // Recieving incremented 2016 from the north:
//        int recieved;
//        MPI_Recv(&recieved, 1, MPI_INT, p_north, 0, cart_comm, MPI_STATUS_IGNORE);
//        printf("Recieving incremented %d from the north: %d\n",sendt, recieved);
//
//        // Sending 2018 westwards
//        sendt = 2018;
//        MPI_Send(&sendt, 1, MPI_INT, p_east, 0, cart_comm);
//    }
//    else if (rank == 3) {
//        int recieved;
//        MPI_Recv(&recieved, 1, MPI_INT, p_west, 0, cart_comm, MPI_STATUS_IGNORE);
//        printf("Recieved from the west: %d\n", recieved);
//    }

    int tag = 0;
    cell *my_test_cell = malloc(10 * sizeof(cell));
    for (int ii = 0; ii < 10; ii++) {
        my_test_cell[ii].strength = ii;
        my_test_cell[ii].color = rank;
    }

    //    SENDING AN ARRAY OF MPI_DERIVED TYPES:
    if (rank >= sqrt(size)) // These ranks have processes to the north:
    {
       MPI_Send(my_test_cell, 10, mpi_cell_t, p_north, tag, cart_comm);
    }
    else
    {
        cell* recieved = malloc(sizeof(cell) * 10);
        MPI_Recv(recieved, 10, mpi_cell_t, p_south, tag, cart_comm, MPI_STATUS_IGNORE);

        for (int i = 0; i < 10; i++)
        {
            printf("RECIEVED FROM SOUTH: %d  %d\n", recieved[i].color, recieved[i].strength);
        }
    }


//    if (rank >= sqrt(size)) // These ranks have processes to the north:
//    {
//        MPI_Send(&rank, 1, MPI_INT, p_north, tag, cart_comm);
//    } else {
//        int recieved;
//        MPI_Recv(&recieved, 1, MPI_INT, p_south, tag, cart_comm, MPI_STATUS_IGNORE);
//        printf("RECIEVED FROM SOUTH: %d\n", recieved);
//    }


    free_multi_dimen_array(array, grid);
    MPI_Finalize();
}

void create_cell_type() {
    // cell type
    const int nitems = 2;
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Aint offsets[2];

    offsets[0] = offsetof(cell, color);
    offsets[1] = offsetof(cell, strength);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cell_t);
    MPI_Type_commit(&mpi_cell_t);
}