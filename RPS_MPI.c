#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

#include "RPS_MPI.h"

void initialize();

void initialize_petri();

void iterate_CA();

void gather_petri();

void create_types();

cell *get_first_row(cell **petri);
cell *get_first_col(cell **petri, int xSize, int ySize);
cell *get_last_row(cell **petri, int ySize);
cell *get_last_col(cell **petri, int xSize, int ySize);
cell** create_full_petri(cell* whole_petri);
cell *get_first_row(cell **petri);
cell *get_first_col(cell **petri, int xSize, int ySize);
cell *get_last_row(cell **petri, int ySize);
cell *get_last_col(cell **petri, int xSize, int ySize);
void print_cell_list(cell *list, int len);
void stitch_left_column(cell *column, int len);
void stitch_right_column(cell *column, int len);
void stitch_top_row(cell *row, int len);
void stitch_bottom_row(cell *row, int len);

void exchange_borders(cell **matrix, int xSize, int ySize, int rank, int size);

int rank;
int size;

// I denote mpi process specific values with hungarian notation, adding a p

// The dimensions of the processor grid. Same for every process
int p_x_dims;
int p_y_dims;

// The location of a process in the process grid. Unique for every process
int p_my_x_dim;
int p_my_y_dim;
int p_local_petri_x_dim;
int p_local_petri_y_dim;

int p_north, p_south, p_east, p_west;

// The dimensions for the process local petri
int p_local_petri_x_dim;
int p_local_petri_y_dim;

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
cell **local_petri_A;
cell **local_petri_B;
cell **petri;

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
    exchange_borders(local_petri_A, p_local_petri_x_dim, p_local_petri_y_dim, rank, size);
    iterate_CA();

    gather_petri();

    // make_bmp(petri, 0);

    MPI_Finalize();

    if (rank == 0) {
        // TODO: Write the petri to an image
    }

    // You should probably make sure to free your memory here
    // We will dock points for memory leaks, don't let your hard work go to waste!
    // free_stuff()

//    for (int i = 0; i < p_local_petri_y_dim; i++) {
//        free(&local_petri_A[i]);
//        free(&local_petri_B[i]);
//    }
//    free(&local_petri_A);
//    free(&local_petri_B);

    exit(0);
}


void create_types() {

    ////////////////////////////////
    ////////////////////////////////
    // cell type
    const int nitems = 2;
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Aint offsets[2];

    offsets[0] = offsetof(cell, color);
    offsets[1] = offsetof(cell, strength);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cell_t);
    MPI_Type_commit(&mpi_cell_t);
    ////////////////////////////////
    ////////////////////////////////



    ////////////////////////////////
    ////////////////////////////////
    // A message for a local petri-dish
    MPI_Type_contiguous(p_local_petri_x_dim * p_local_petri_y_dim,
                        mpi_cell_t,
                        &local_petri_t);
    MPI_Type_commit(&local_petri_t);
    ////////////////////////////////
    ////////////////////////////////


    //TODO: Create MPI types for border exchange

}


void initialize() {
    // The dimensions for the process local petri
    int square = sqrt(size);
    p_local_petri_x_dim = (IMG_X / square) + 2; // +2 for the borders on each side.
    p_local_petri_y_dim = (IMG_Y / square) + 2; // +2 for the borders on each side.

    // TODO: When allocating these buffers, keep in mind that you might need to allocate a little more
    // than just your piece of the petri.
    local_petri_A = malloc(p_local_petri_y_dim * sizeof(cell *));
    local_petri_B = malloc(p_local_petri_y_dim * sizeof(cell *));

    for (int i = 0; i < p_local_petri_x_dim; i++) {
        local_petri_A[i] = malloc(p_local_petri_x_dim * sizeof(cell));
        local_petri_B[i] = malloc(p_local_petri_x_dim * sizeof(cell));
    }

    int xstart = 0;
    int ystart = 0;

    // TODO: Randomly the local dish. Only perturb cells that belong to your process,
    // Seed some CAs
    for (int ii = 0; ii < 100 / size; ii++) {
        int rx = rand() % (p_local_petri_x_dim - 2);
        int ry = rand() % (p_local_petri_y_dim - 2);
        int rt = rand() % 4;

        local_petri_A[rx][ry].color = rt;
        local_petri_A[rx][ry].strength = 1;
    }
}

void iterate_CA() {
    iterate_image2(local_petri_A, local_petri_B, p_local_petri_x_dim, p_local_petri_y_dim);
}

void gather_petri() {
    int grid = (p_local_petri_x_dim * p_local_petri_y_dim);
    cell *petri_package = malloc(sizeof(cell) * grid);
    int i = 0;
    for (int y = 0; y < p_local_petri_y_dim; y++) {
        for (int x = 0; x < p_local_petri_x_dim; x++) {
            petri_package[i++] = local_petri_B[y][x];
        }
    }

    int sq = sqrt(size);
    if (rank == 0) {
        cell* whole_petri = malloc(sizeof(cell) * ((IMG_Y + (sq*2))*( IMG_X + (sq*2))));
        MPI_Gather(&petri_package, grid, mpi_cell_t, whole_petri, grid, mpi_cell_t, 0, cart_comm);
        create_full_petri(whole_petri);
    } else {
        MPI_Gather(&petri_package, grid, mpi_cell_t, NULL, grid, mpi_cell_t, 0, cart_comm);
    }

}

cell** create_full_petri(cell* whole_petri)
{
    // GRID SIZE TO REDUCE CODE LENGTH.
    int grid = ((p_local_petri_x_dim/size) * (p_local_petri_y_dim/size));
    int sq = sqrt(size);
    int empty_cell_offset = (sq*2);

    // ALLOCATING SPACE FOR A 2D ARRAY OF THE ENTIRE PETRI:
    petri = malloc(sizeof(*petri)*IMG_Y);
    for (int y = 0; y < IMG_Y; y++)
        petri[y] = malloc(sizeof(*petri)*IMG_X);
//    printf("size of whole grid = %d\n", (IMG_Y) * (IMG_X ));


//    printf("Local x dimension: %d\n", p_local_petri_x_dim);
//    printf("Local y dimension: %d\n", p_local_petri_y_dim);

    int i = 2;
    for (int s = size-1; s >= 0 ; s--) {
        int mod = ((size -1) - s);

        // Setting up the next iteration over the big array:
        for (int y = 1; y < p_local_petri_y_dim-1; y++) {
            for (int x = 1; x < p_local_petri_x_dim-1; x++) {
                petri[y][x] = whole_petri[i];
                printf("MAPPING: (%d, %d) == %d    CONTENT: [%d, %d]\n", x, y, i, petri[x][y].color, petri[x][y].strength);
                i++;
            }
        }
        i++;
    }
}

cell *get_first_row(cell **petri) {
    return petri[1];
}

cell *get_last_row(cell **petri, int ySize) {
    return petri[ySize - 2];
}

cell *get_first_col(cell **petri, int xSize, int ySize) {
    cell *output = malloc(ySize * sizeof(cell));
    for (int y = 0; y < ySize; y++) {
        output[y] = petri[y][1];
    }
    return output;
}

cell *get_last_col(cell **petri, int xSize, int ySize) {
    cell *output = malloc(ySize * sizeof(cell));
    for (int y = 0; y < ySize; y++) {
        output[y] = petri[y][xSize - 2];
    }
    return output;
}

void exchange_borders(cell **matrix, int xSize, int ySize, int rank, int size) {
    int squared = sqrt(size);

    if (0 <= rank + 1 && rank + 1 <= squared) {
//        printf("Rank %d is first row, exchange downwards.\n", rank);

        // SEND PACKAGE SOUTH ONLY:
        cell *south_package = get_first_row(matrix);
        MPI_Send(&south_package, xSize, mpi_cell_t, p_south, 0, cart_comm);

        // RECIEVE PACKAGE FROM SOUTH ONLY:
        cell *recieve_south = malloc(xSize * sizeof(cell));
        MPI_Recv(recieve_south, xSize, mpi_cell_t, p_south, 0, cart_comm, MPI_STATUS_IGNORE);
        stitch_bottom_row(recieve_south, xSize);

    } else if ((squared * (squared - 1)) + 1 <= rank + 1 && rank + 1 <= size) {
//        printf("Rank %d is last row. exchange upwards.\n", rank);
        // SEND PACKAGE NORTH ONLY:
        cell *north_package = get_last_row(matrix, ySize);
        MPI_Send(&north_package, xSize, mpi_cell_t, p_north, 0, cart_comm);

        // RECIEVE PACKAGE FROM NORTH ONLY:
        cell *recieve_north = malloc(xSize * sizeof(cell));
        MPI_Recv(recieve_north, xSize, mpi_cell_t, p_north, 0, cart_comm, MPI_STATUS_IGNORE);
        stitch_top_row(recieve_north, xSize);

    } else {
        // SEND PACKAGE BOTH NORTH AND SOUTH ONLY:
        cell *south_package = get_first_row(matrix);
        MPI_Send(&south_package, xSize, mpi_cell_t, p_south, 0, cart_comm);

        cell *north_package = get_last_row(matrix, ySize);
        MPI_Send(&north_package, xSize, mpi_cell_t, p_north, 0, cart_comm);

        // RECIEVE PACKAGE BOTH NORTH AND SOUTH ONLY:
        cell *recieve_south = malloc(xSize * sizeof(cell));
        MPI_Recv(recieve_south, xSize, mpi_cell_t, p_south, 0, cart_comm, MPI_STATUS_IGNORE);
        stitch_bottom_row(recieve_south, xSize);

        cell *recieve_north = malloc(xSize * sizeof(cell));
        MPI_Recv(recieve_north, xSize, mpi_cell_t, p_north, 0, cart_comm, MPI_STATUS_IGNORE);
        stitch_top_row(recieve_north, xSize);
    }

    bool notFound = true;
    for (int i = 0; i < squared; i++) {
        if (rank == squared * i) {
//            printf("Rank %d is first column, exchange left.\n", rank);
            notFound = false;
            cell *west_package = get_first_col(matrix, ySize, ySize);
            MPI_Send(&west_package, ySize, mpi_cell_t, p_west, 0, cart_comm);

            // RECIEVE PACKAGE BOTH NORTH AND SOUTH ONLY:
            cell *recieve_east = malloc(ySize * sizeof(cell));
            MPI_Recv(recieve_east, ySize, mpi_cell_t, p_west, 0, cart_comm, MPI_STATUS_IGNORE);
            stitch_right_column(recieve_east, ySize);
            break;
        } else if (rank == (squared * (i + 1)) - 1) {
//            printf("Rank %d is last column. exchange right.\n", rank);
            notFound = false;
            cell *east_package = get_last_col(matrix, xSize, ySize);
            MPI_Send(&east_package, ySize, mpi_cell_t, p_east, 0, cart_comm);

            // RECIEVE PACKAGE BOTH NORTH AND SOUTH ONLY:
            cell *recieve_west = malloc(ySize * sizeof(cell));
            MPI_Recv(recieve_west, ySize, mpi_cell_t, p_east, 0, cart_comm, MPI_STATUS_IGNORE);
            stitch_left_column(recieve_west, ySize);
            break;
        }
    }
    if (notFound) {
//        printf("Rank %d is middle column, exchange right and left.\n", rank);

        // RECIEVE PACKAGE BOTH NORTH AND SOUTH ONLY:
        cell *west_package = get_first_col(matrix, ySize, ySize);
        MPI_Send(&west_package, ySize, mpi_cell_t, p_west, 0, cart_comm);
        cell *east_package = get_last_col(matrix, xSize, ySize);
        MPI_Send(&east_package, ySize, mpi_cell_t, p_east, 0, cart_comm);

        // RECIEVE PACKAGE BOTH NORTH AND SOUTH ONLY:
        cell *recieve_east = malloc(ySize * sizeof(cell));
        MPI_Recv(recieve_east, ySize, mpi_cell_t, p_west, 0, cart_comm, MPI_STATUS_IGNORE);
        cell *recieve_west = malloc(ySize * sizeof(cell));
        MPI_Recv(recieve_west, ySize, mpi_cell_t, p_east, 0, cart_comm, MPI_STATUS_IGNORE);

        stitch_right_column(recieve_east, ySize);
        stitch_left_column(recieve_west, ySize);
    }
}

void stitch_bottom_row(cell *row, int len) {
    for (int i = 0; i < len; i++) {
        local_petri_A[p_local_petri_y_dim - 1][i] = row[i];
    }
}

void stitch_top_row(cell *row, int len) {
    for (int i = 0; i < len; i++) {
        local_petri_A[0][i] = row[i];
    }
}

void stitch_left_column(cell *column, int len) {
    for (int i = 0; i < len; i++) {
        local_petri_A[i][0] = column[i];
    }
}

void stitch_right_column(cell *column, int len) {
    for (int i = 0; i < len; i++) {
        local_petri_A[i][p_local_petri_x_dim - 1] = column[i];
    }
}

void print_cell_list(cell *list, int len) {
    for (int i = 0; i < len; i++) {
        printf("{%d, %d} ", list[i].color, list[i].strength);
    }
    printf("\n");
}
