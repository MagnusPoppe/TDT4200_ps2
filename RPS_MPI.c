#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

#include "RPS_MPI.h"

int rank;
int size;

// The dimensions of the processor grid. Same for every process
int p_x_dims;
int p_y_dims;

// The location of a process in the process grid. Unique for every process
int p_my_x_dim;
int p_my_y_dim;

int p_north, p_south, p_east, p_west;

// The dimensions for the process local petri
int local_x;
int local_y;

MPI_Comm cart_comm;

// some datatypes, useful for sending data with somewhat less primitive semantics
MPI_Datatype border_row_t;  // TODO: Implement this
MPI_Datatype border_col_t;  // TODO: Implement this
MPI_Datatype local_petri_t; // Already implemented
MPI_Datatype mpi_cell_t;    // Already implemented

// Each process is responsible for one part of the petri dish.
// Since we can't update the petri-dish in place each process actually
// gets two petri-dishes which they update in a lockstep fashion.
// dish A is updated by writing to dish B, then next step dish B updates dish A.
// (or you can just swap them inbetween iterations)
int **local_petri_A;
int **local_petri_B;
int **petri;


int main(int argc, char **argv) {

    srand(1234);

    // Ask MPI what size (number of processors) and rank (which process we are)
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    ////////////////////////////////
    // Create cartesian communicator
    int dims[2];
    dims[0] = p_x_dims;
    dims[1] = p_y_dims;

    int periods[2]; // we set these to 0 because we are not interested in wrap-around
    periods[0] = 0;
    periods[1] = 0;

    int coords[2];
    coords[0] = p_my_x_dim;
    coords[1] = p_my_y_dim;

    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    MPI_Cart_shift(cart_comm, 0, 1, &p_north, &p_south);
    MPI_Cart_shift(cart_comm, 1, 1, &p_west, &p_east);

    p_x_dims = dims[0];
    p_y_dims = dims[1];

    p_my_x_dim = coords[0];
    p_my_y_dim = coords[1];


    // RUNNING ALGORITHM:
    initialize();
    create_types();
    exchange_borders(
            local_petri_A, local_x, local_y, rank, size,
            p_north, p_south, p_east, p_west, border_row_t, border_col_t, cart_comm
    );
//    iterate_CA();
    gather_petri();

    // make_bmp(petri, 0);

    MPI_Finalize();

    if (rank == 0) {
        // TODO: Write the petri to an image
    }

    // You should probably make sure to free your memory here
    // We will dock points for memory leaks, don't let your hard work go to waste!
    // free_stuff()

//    for (int i = 0; i < local_y; i++) {
//        free(&local_petri_A[i]);
//        free(&local_petri_B[i]);
//    }
//    free(&local_petri_A);
//    free(&local_petri_B);

    exit(0);
}

void create_types() {

    // cell type
    const int nitems = 2;
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Aint offsets[2];

    offsets[0] = offsetof(cell, color);
    offsets[1] = offsetof(cell, strength);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cell_t);
    MPI_Type_commit(&mpi_cell_t);

    // A message for a local petri-dish
    MPI_Type_contiguous(local_x * local_y,
                        MPI_INT,
                        &local_petri_t);
    MPI_Type_commit(&local_petri_t);

    // MESSAGES FOR BORDER EXCHANGE
    MPI_Type_contiguous(local_x,
                        MPI_INT,
                        &border_row_t);
    MPI_Type_commit(&border_row_t);

    MPI_Type_contiguous(local_y,
                        MPI_INT,
                        &border_col_t);
    MPI_Type_commit(&border_col_t);
}

void initialize() {
    // The dimensions for the process local petri
    int square = sqrt(size);
    local_x = (IMG_X / square) + 2; // +2 for the borders on each side.
    local_y = (IMG_Y / square) + 2; // +2 for the borders on each side.

    // TODO: When allocating these buffers, keep in mind that you might need to allocate a little more
    // than just your piece of the petri.
    local_petri_A = malloc(local_y * sizeof(int *));
    local_petri_B = malloc(local_y * sizeof(int *));

    int c = 0;
    if (rank == 0) printf("SQUARE + BORDERS OF THE GRID. BORDER SIZE == 1\n");
    for (int i = 0; i < local_y; i++) {
        local_petri_A[i] = malloc(local_x * sizeof(int));
        local_petri_B[i] = malloc(local_x * sizeof(int));
        for (int x = 0; x < local_x; x++) {
        local_petri_A[i][x] = c;
        local_petri_B[i][x] = c++;
        }
        if (rank == 0) print_int_list(local_petri_A[i], local_x);

    }

//    int xstart = 0;
//    int ystart = 0;
//
//    // TODO: Randomly the local dish. Only perturb ints that belong to your process,
//    // Seed some CAs
//    for (int ii = 0; ii < 100 / size; ii++) {
//        int rx = rand() % (local_x - 2);
//        int ry = rand() % (local_y - 2);
//        int rt = rand() % 4;
//
//        local_petri_A[rx][ry].color = rt;
//        local_petri_A[rx][ry].strength = 1;
//    }
}

//void iterate_CA() {
//    iterate_image2(local_petri_A, local_petri_B, local_x, local_y);
//}

void gather_petri() {
    int grid = (local_x * local_y);
    int *petri_package = malloc(sizeof(int) * grid);
    int i = 0;
    for (int y = 0; y < local_y; y++) {
        for (int x = 0; x < local_x; x++) {
            petri_package[i++] = local_petri_B[y][x];
        }
    }

    int sq = sqrt(size);
    if (rank == 0) {
        int* whole_petri = malloc(sizeof(int) * ((IMG_Y + (sq*2))*( IMG_X + (sq*2))));
        MPI_Gather(&petri_package, grid, MPI_INT, whole_petri, grid, MPI_INT, 0, cart_comm);
        create_full_petri(whole_petri);
    } else {
        MPI_Gather(&petri_package, grid, MPI_INT, NULL, grid, MPI_INT, 0, cart_comm);
    }

}

int** create_full_petri(int* whole_petri)
{
    // GRID SIZE TO REDUCE CODE LENGTH.
    int grid = ((local_x/size) * (local_y/size));
    int sq = sqrt(size);
    int empty_int_offset = (sq*2);

    // ALLOCATING SPACE FOR A 2D ARRAY OF THE ENTIRE PETRI:
    petri = malloc(sizeof(*petri)*IMG_Y);
    for (int y = 0; y < IMG_Y; y++)
        petri[y] = malloc(sizeof(*petri)*IMG_X);
//    printf("size of whole grid = %d\n", (IMG_Y) * (IMG_X ));


//    printf("Local x dimension: %d\n", local_x);
//    printf("Local y dimension: %d\n", local_y);

    int i = 2;
    for (int s = size-1; s >= 0 ; s--) {
        int mod = ((size -1) - s);

        // Setting up the next iteration over the big array:
        for (int y = 1; y < local_y-1; y++) {
            for (int x = 1; x < local_x-1; x++) {
                petri[y][x] = whole_petri[i];
//                printf("MAPPING: (%d, %d) == %d    CONTENT: [%d, %d]\n", x, y, i, petri[x][y].color, petri[x][y].strength);
                i++;
            }
        }
        i++;
    }
}

void print_int_list(int *list, int len) {
    for (int i = 0; i < len; i++) {
        if (list[i] < 10) printf("0%d ", list[i]);
        else              printf("%d ",  list[i]);
    }
    printf("\n");
}
