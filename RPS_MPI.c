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
int EXTRA_EDGE = 4;

MPI_Comm cart_comm;

// some datatypes, useful for sending data with somewhat less primitive semantics
MPI_Datatype border_row_t;
MPI_Datatype border_col_t;
MPI_Datatype local_petri_t;
MPI_Datatype mpi_cell_t;

// Each process is responsible for one part of the petri dish.
// Since we can't update the petri-dish in place each process actually
// gets two petri-dishes which they update in a lockstep fashion.
// dish A is updated by writing to dish B, then next step dish B updates dish A.
// (or you can just swap them inbetween iterations)
cell **local_petri_A;
cell **local_petri_B;
cell **petri;
cell*** images;


int main(int argc, char **argv) {


    // Ask MPI what size (number of processors) and rank (which process we are)
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(1234 * rank);

    // RUNNING ALGORITHM:
    create_cartesian_communicator();
    initialize();
    create_types();

    if (rank == 0)
    {
        // THIS IS DEALLOCATED INSIDE THE "WRITE IMAGES" LOOP.
        images = calloc(ITERATIONS, sizeof(cell**));
    }

    for( int i = 0; i < ITERATIONS; i++) {
        exchange_borders(
            local_petri_A, local_x, local_y, rank, size,
            p_north, p_south, p_east, p_west, border_row_t, border_col_t, cart_comm
        );

        iterate_CA(local_petri_A, local_petri_B);
        gather_petri();
        if (rank == 0)
        {
            images[i] = petri;
        }

        // FREEING THE LAST USED PETRI DISH.
        free_multi_cell_array(local_petri_A, local_x);
        local_petri_A = local_petri_B;

        local_petri_B = calloc(local_y, sizeof(cell *));
        for (int i = 0; i < local_y; i++)
            local_petri_B[i] = calloc(local_x, sizeof(cell));
    }
    // SINCE THE POINTERS ALWAYS CHANGE, ALL LOCAL PETRI B WILL BE DE-ALLOCATED AS LOCAL PETRI A.
    // DE ALLOCATING THE LAST ONE.
    free_multi_cell_array(local_petri_B, local_x);

    // WRITING IMAGE AND FREEING THE IMAGES FROM MEMORY. FREEING THESE IMAGES ALSO FREES
    // ALL "PETRI" MULTI DIMENSIONAL ARRAYS ALLOCATED EARLIER IN THE APP.
    write_images();
    MPI_Finalize();

    // You should probably make sure to free your memory here
    // We will dock points for memory leaks, don't let your hard work go to waste!
    // free_stuff()

    exit(0);
}

void write_images()
{
    if (rank == 0) {
        for (int i = 0; i < ITERATIONS; i++) {
            make_bmp(images[i], i);
            free_multi_cell_array(images[i], IMG_X);
        }
    }
    free(images);
}

void create_cartesian_communicator() {

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
                        mpi_cell_t,
                        &local_petri_t);
    MPI_Type_commit(&local_petri_t);

    // MESSAGES FOR BORDER EXCHANGE
    MPI_Type_contiguous(local_x,
                        mpi_cell_t,
                        &border_row_t);
    MPI_Type_commit(&border_row_t);

    MPI_Type_contiguous(local_y,
                        mpi_cell_t,
                        &border_col_t);
    MPI_Type_commit(&border_col_t);
}

void initialize() {
    // The dimensions for the process local petri
    int square = sqrt(size);
    local_x = (IMG_X / square) + EXTRA_EDGE;
    local_y = (IMG_Y / square) + EXTRA_EDGE;

    // TODO: When allocating these buffers, keep in mind that you might need to allocate a little more
    // than just your piece of the petri.
    local_petri_A = calloc(local_y, sizeof(cell *));
    local_petri_B = calloc(local_y, sizeof(cell *));

    for (int i = 0; i < local_y; i++) {
        local_petri_A[i] = calloc(local_x, sizeof(cell));
        local_petri_B[i] = calloc(local_x, sizeof(cell));
    }

    // Seed some CAs
    for (int ii = 0; ii < 100; ii++) {
        int rx = rand() % (local_x);
        int ry = rand() % (local_y);
        int rt = rand() % 4;

        local_petri_A[rx][ry].color = rt;
        local_petri_A[rx][ry].strength = 1;
    }
}

void iterate_CA(cell** current_image, cell** next_image) {
    iterate_image2(current_image, next_image, local_x, local_y);
}

void gather_petri() {
    int grid = (local_x * local_y);

    // Allocating memory for payload. Freed at the end of the function.
    cell *petri_package = malloc(sizeof(cell) * grid);
    int i = 0;
    for (int y = 0; y < local_y; y++) {
        for (int x = 0; x < local_x; x++) {
            petri_package[i++] = local_petri_B[y][x];
        }
    }

    // Sending and recieving using MPI_Gather()
    int sq = sqrt(size);
    if (rank == 0)
    {
        // Recieving data. Allocating memory for the whole image. Deallocated at the end of IF.
        cell* whole_petri = malloc(sizeof(cell) * ((IMG_Y + (sq*EXTRA_EDGE))*( IMG_X + (sq*EXTRA_EDGE))));
        MPI_Gather(petri_package, 1, local_petri_t, whole_petri, 1, local_petri_t, 0, MPI_COMM_WORLD);

        // Formatting recieved data:
        convert_1D_to_2D(whole_petri, ((IMG_Y + (sq*EXTRA_EDGE))*( IMG_X + (sq*EXTRA_EDGE))));

        // Freeing recieved data:
        free(whole_petri);
    }
    else {
        MPI_Gather(petri_package, 1, local_petri_t, NULL, 1, local_petri_t, 0, MPI_COMM_WORLD);
    }
    // Freeing payload
    free(petri_package);
}

void convert_1D_to_2D(cell *whole_petri, int len) {
    int sq = sqrt(size);

    // ALLOCATING MEMORY FOR ARRAY TO SORT ALL CELLS. FREED AT THE END OF LOOP.
    cell** temppetri = malloc(sizeof(cell*) * (IMG_Y + (sq*EXTRA_EDGE)));
    for (int i = 0; i < IMG_X + sq*EXTRA_EDGE; i++)
        temppetri[i] = malloc(sizeof(cell) * ( IMG_X + (sq*EXTRA_EDGE)) ) ;

    // ALLOCATING MEMORY FOR MAIN IMAGE. THIS IS DEALLOCATED AFTER IMAGES ARE SAVED.
    petri = malloc(sizeof(cell*) * IMG_Y);
    for (int i = 0; i < IMG_X; i++)
        petri[i] = malloc(sizeof(cell) * IMG_X);


    int index = 0;
    for (int i = 0; i < size; i++) {
        int xstart = (i % sq) * local_x;
        int ystart = (i / sq) * local_x;
        for (int y = ystart; y < ystart + local_y; y++) {
            for (int x = xstart; x < xstart + local_x; x++) {
                temppetri[y][x] = whole_petri[index];
                index++;
            }
        }
    }

    int yreal = 0;
    int ystride = 0;
    for (int y = EXTRA_EDGE/2; y < ( IMG_Y + (sq*EXTRA_EDGE)) - (EXTRA_EDGE/2); y++)
    {
        int xreal = 0;
        int xstride = 0;
        for (int x = EXTRA_EDGE/2; x < ( IMG_X + (sq*EXTRA_EDGE)) - (EXTRA_EDGE/2); x++)
        {
            if (xstride <= IMG_X / sq -1 && ystride <= IMG_Y / sq -1)  {
                petri[yreal][xreal] = temppetri[y][x];
                xreal++;
                xstride++;
            } else {
                x += EXTRA_EDGE/2 +1;
                xstride = 0;
            }
        }
        if (ystride <= IMG_Y / sq -1)  {
            ystride++;
            yreal++;
        } else {
            y += EXTRA_EDGE/2 +1;
            ystride = 0;
        }
    }
    free_multi_cell_array(temppetri, ( IMG_X + (sq*EXTRA_EDGE)));
}

void free_multi_cell_array(cell **array, int width)
{
    for (int i = 0; i < width; i++)
        free(array[i]);
    free(array);
}
