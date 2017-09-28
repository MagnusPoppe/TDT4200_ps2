//
// Created by Magnus Poppe Wang on 27/09/2017.
//

#include "border_exchange.h"

cell* get_first_row(cell** petri)
{
    return petri[1];
}
cell* get_last_row(cell** petri, int ySize)
{
    return petri[ySize-2];
}
cell* get_first_col(cell **petri, int xSize, int ySize)
{
    return NULL;
}
cell* get_last_col(cell **petri,  int xSize, int ySize)
{
    return NULL;
}


void exchange_borders(cell** matrix, int xSize, int ySize, int rank, int size)
{
    //TODO: Exchange borders inbetween each step
    int squared = sqrt(size);

    if ( 0 <= rank+1 && rank+1 <= squared )  // Øverste rad
    {
        printf("Rank %d is first row, exchange downwards.\n", rank);

        // SEND PACKAGE SOUTH ONLY:
        cell* south_package = get_first_row(matrix);
        MPI_Send(&south_package, xSize, mpi_cell_t, p_south, 0, cart_comm);
        free(south_package);

        // RECIEVE PACKAGE FROM SOUTH ONLY:
        cell* recieve_south = malloc(xSize* sizeof(cell));
        MPI_Recv(&recieve_south, xSize, mpi_cell_t, p_south, 0, cart_comm, MPI_STATUS_IGNORE);

        for (int i = 0; i < xSize; i++)
        {
            printf("[%d, %d]\n", recieve_south[i].strength, recieve_south[i].color);
        }
    }
    else if ( (squared*(squared-1)) +1 <= rank+1 && rank+1 <= size )  // Nederste rad
    {
        printf("Rank %d is last row. exchange upwards.\n", rank);
        // SEND PACKAGE NORTH ONLY:
        cell* north_package = get_last_row(matrix, ySize);
        MPI_Send(&north_package, xSize, mpi_cell_t, p_north, 0, cart_comm);
        free(north_package);

        // RECIEVE PACKAGE FROM NORTH ONLY:
        cell* recieve_north = malloc(xSize* sizeof(cell));
        MPI_Recv(&recieve_north, xSize, mpi_cell_t, p_north, 0, cart_comm, MPI_STATUS_IGNORE);

        for (int i = 0; i < xSize; i++)
        {
            printf("[%d, %d]\n", recieve_north[i].strength, recieve_north[i].color);
        }
    }
//    else
//    {
//        cell* south_package = get_first_row(local_petri_A);
//        cell* north_package = get_last_row(local_petri_A, p_local_petri_y_dim);
//
//        // SEND ONE PACKAGE NORTH AND ONE SOUTH:
//        MPI_Send(&north_package, p_local_petri_x_dim, mpi_cell_t, p_north, 0, cart_comm);
//        MPI_Send(&south_package, p_local_petri_x_dim, mpi_cell_t, p_south, 0, cart_comm);
//    }

    bool notFound = true;
    for (int i = 0; i < squared; i++)
    {
        if ( rank == squared*i) {
            printf("Rank %d is first column, exchange right.\n", rank);
            notFound = false;
            break;
        }
        else if ( rank == (squared*(i+1)) -1 )
        {
            printf("Rank %d is last column. exchange left.\n", rank);
            notFound = false;
            break;
        }
    }
    if (notFound) printf("Rank %d is middle column, exchange right and left.\n", rank);
}