//
// Created by Magnus Poppe Wang on 27/09/2017.
//

#include "help.h"
void free_multi_dimen_array(int** array, int size)
{
    for (int i = 0; i < size; i++) {
        free(array[i]);
    }
    free(array);
}

int** create_2D_array(int size) {
    // Allocating array space in memory:
    int** array = malloc(size* sizeof(int*));
    for (int i = 0; i < size; i++) {
        array[i] = malloc(size * sizeof(int));
    }

    // Filling array with ints;
    int i = 0;
    for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
            array[y][x] = i++;
        }
    }
    return array;
}

void print_array(int** array, int size) {
    for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
            if (array[y][x] < 10) printf("0%d ", array[y][x]);
            else printf("%d ", array[y][x]);
        }
        printf("\n");
    }
}