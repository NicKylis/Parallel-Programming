#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>
#include <cblas.h>
#include "knn_cilk_pthreads.h"
#include <time.h>
//#include <hdf5.h>
#include <pthread.h>
#include <cilk/cilk.h>

/*
    !!WARNING!! Before testing this code I suggest testing the omp code, there are more detailed comments there

    This code performs a knn search as instructed. I have created multiple knn functions as you will see
    in the header file, the fastest one being knnsearch_parallel and the most simple (and accurate) one being knnsearch and knnsearchBLAS.
    More details in the header file.
    This implementation uses a .txt file to read data resulting in a loss in accuracy. At the bottom of this page,
    I proivde the experimental code that couldn't work trying to read the data of an hdf5 file as instructed.
    As a final note, this code works best for C=Q, but also works for C!=Q.
*/

// To compile you will need (default): /yourPath/clang -I/usr/include/openblas -L/usr/lib/openblas -o cilk_pthreads cilk_pthreads.c -lopenblas -lm -fopencilk -O3

int main() {
    // srand(100); // For random value initialization
    // clock_t s1, s2, s3; // Uncomment these if you need timespamps 
    // double cpu1, cpu2;

    // Parameters
    int num_points = 10000;      // Number of points in C and Q
    int dimensions = 32;         // Dimensionality of each point
    int k = 10;                  // Number of nearest neighbors to find

    // Allocate memory for C and Q
    double* C = (double*)malloc(num_points * dimensions * sizeof(double));
    double* Q = (double*)malloc(num_points * dimensions * sizeof(double));

    FILE *file = fopen("arrays.txt", "r"); // Write the filename.txt here (default is provided)
    // Note: if another file is used, it should be noted that this implementation ignores the header line of the .txt
    // The following code reads two arrays so if there is a segmantation fault it might be that the .txt file had one
    
    char header[100];
    fgets(header, sizeof(header), file);

     for (int i = 0; i < num_points * dimensions; i++) {
        if (fscanf(file, "%lf %lf", &C[i], &Q[i]) != 2) {
            fprintf(stderr, "Error reading data at line %d\n", i + 2);  // +2 to account for the header line
            fclose(file);
            return EXIT_FAILURE;
        }
    }

    // Close the file
    fclose(file);

    // Example initialization (of random data)
    // for (int i = 0; i < num_points * dimensions; ++i) {
    //     C[i] = rand() / (double)RAND_MAX;
    //     Q[i] = rand() / (double)RAND_MAX;
    // }


    // Allocate memory for output arrays
    int* indices = (int*)malloc(num_points * k * sizeof(int)); // Output indices
    double* distances = (double*)malloc(num_points * k * sizeof(double)); // Output distances

    // Call the knnsearch function (I recommend to at least check out knnsearch and knnsearch_parallel)
    // Note: on all implementations of the knnsearch, any C,Q input combination work - with different accuracy
    knnsearch(C, C, num_points, dimensions, k, indices, distances);

    // Example output (print the first point's k-nearest neighbors)
    printf("Nearest neighbors for the first query point:\n");
    for (int i = 0; i < k; ++i) {
        printf("Index: %d, Distance: %f\n", indices[i], distances[i]);
    }

    // Free allocated memory
    free(C);
    free(Q);
    free(indices);
    free(distances);

    return 0;
}
