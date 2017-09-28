//
// Created by Magnus Poppe Wang on 24/09/2017.
//


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// #include "RPS_MPI.h"


// void experiment_multi_dimens_array();
// void test_local_petri_dish();
int** test_multi_dimens_array();
void print_array(int** array);
void test_exchange_borders(int size);
void exchange_borders(int rank, int size);

int main(int argc, char** argv) {
//    // experiment_multi_dimens_array();
//    int size = 10;
//
//    // for( int rank = 0; rank < size; rank++ ) {
//    //     test_local_petri_dish(rank, size);
//    //     printf("\n\n");
//    // }
//
//    int** array = create_2D_array(size);
//
//    print_array(array);
//
//    free_multi_dimen_array(array, size);

    // TEST BORDERING EDGES OF AN ARRAY:

    // test_exchange_borders(4);
}

void test_exchange_borders(int size)
{
    printf("Size         = %d\n", size);
    printf("Squared size = %d\n", (int)(sqrt(size)));
    for (int rank = 0; rank < size; rank++)
    {
        exchange_borders(rank, size);
    }
}

void exchange_borders(int rank, int size)
{
    //TODO: Exchange borders inbetween each step
    int squared = sqrt(size);

    if ( 0 <= rank+1 && rank+1 <= squared )  // Øverste rad
    {
        printf("Rank %d is first row, exchange downwards.\n", rank);
    }
    else if ( (squared*(squared-1)) +1 <= rank+1 && rank+1 <= size )  // Nederste rad
    {
        printf("Rank %d is last row. exchange upwards.\n", rank);
    }
    else
    {
        printf("Rank %d is middle row, exchange up and down.\n", rank);
    }

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

//void free_multi_dimen_array(int** array, int size)
//{
//    for (int i = 0; i < size; i++) {
//        free(array[i]);
//    }
//    free(array);
//}
//
//int** create_2D_array(int size) {
//    int** array = malloc(size* sizeof(int*));
//
//    for (int i = 0; i < size; i++) {
//        array[i] = malloc(size* sizeof(int));
//    }
//
//    int i = 0;
//    for (int y = 0; y < 10; y++) {
//        for (int x = 0; x < 10; x++) {
//            array[y][x] = i++;
//        }
//    }
//    return array;
//}
//
//void print_array(int** array) {
//    for (int y = 0; y < 10; y++) {
//        for (int x = 0; x < 10; x++) {
//            if (array[y][x] < 10) printf("0%d ", array[y][x]);
//            else printf("%d ", array[y][x]);
//        }
//        printf("\n");
//    }
//}
//
//void test_local_petri_dish(int rank, int size)
//{
//
//
//    // The dimensions for the process local petri
//    int p_local_petri_x_dim = IMG_X;
//    int p_local_petri_y_dim;
//
//    //TODO: assign the following to something more useful than 0
//    p_local_petri_y_dim = (IMG_Y / size); // +2 for the borders on each side.
//
//    // TODO: When allocating these buffers, keep in mind that you might need to allocate a little more
//    // than just your piece of the petri.
//    cell** local_petri_A = malloc(p_local_petri_y_dim*sizeof(cell*));
//    cell** local_petri_B = malloc(p_local_petri_y_dim*sizeof(cell*));
//
//    for (int i = 0; i < p_local_petri_x_dim; i++)
//    {
//        local_petri_A[i] = malloc(p_local_petri_x_dim * sizeof(cell));
//        local_petri_B[i] = malloc(p_local_petri_x_dim * sizeof(cell));
//    }
//
//    int xstart = 0;
//    int ystart = p_local_petri_y_dim * rank;
//
//    printf("Local petri dish now created.\n");
//    printf("Domain of rank %d:\n", rank);
//    printf("x1=%d \tx2=%d:\n", xstart, (xstart+p_local_petri_x_dim));
//    printf("y1=%d \ty2=%d:\n", ystart, (ystart+p_local_petri_y_dim));
//
//    for (int i = 0; i < p_local_petri_y_dim; i++)
//    {
//        free(local_petri_A[i]);
//        free(local_petri_B[i]);
//    }
//    // TODO: Randomly the local dish. Only perturb cells that belong to your process,
//    // Seed some CAs
//    for(int ii = 0; ii < 100/size; ii++){
//        int rx = rand() % (p_local_petri_x_dim - 1);
//        int ry = rand() % (p_local_petri_y_dim - 1);
//        int rt = rand() % 4;
//
//        local_petri_A[rx][ry].color = rt;
//        local_petri_A[rx][ry].strength = 1;
//        printf("NEW CELL (%d, %d) with color %d and strength %d\n", rx, ry, rt, 1);
//    }
//}
//
//void experiment_multi_dimens_array()
//{
//    // int* x_size = malloc(5 * sizeof(int));
//
//    int y_size = 10;
//    int x_size = 400;
//
//    int** tot_size = malloc(y_size * sizeof(int*));
//
//    for (int i = 0; i < y_size; i++)
//    {
//        tot_size[i] = malloc(x_size * sizeof(int));
//    }
//
//    for (int y = 0; y < y_size; y++)
//    {
//        for (int x = 0; x < x_size; x++)
//        {
//            tot_size[y][x] = x+y;
//        }
//    }
//
//    for (int y = 0; y < y_size; y++)
//    {
//        for (int x = 0; x < x_size; x++)
//        {
//            printf("%d ", tot_size[y][x]);
//        }
//        printf("\n");
//    }
//
//    for (int i = 0; i < y_size; i++)
//    {
//        free(tot_size[i]);
//    }
//}