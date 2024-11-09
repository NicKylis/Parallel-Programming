#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <cblas.h>
#include "functions.h"
#include <time.h>


int main() {
    srand(100);
    clock_t s1, s2, s3;
    double cpu1, cpu2;
    // Example parameters
    int num_points = 1000;      // Number of points in C and Q
    int dimensions = 32;         // Dimensionality of each point
    int k = 10;                  // Number of nearest neighbors to find

    // Allocate memory for C and Q
    double* C = (double*)malloc(num_points * dimensions * sizeof(double));
    double* Q = (double*)malloc(num_points * dimensions * sizeof(double));

    FILE *file = fopen("two_arrays_32000_points.txt", "r");
    
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

    s1 = clock();
    // Example initialization (of random data)
    // for (int i = 0; i < num_points * dimensions; ++i) {
    //     C[i] = rand() / (double)RAND_MAX;
    //     Q[i] = rand() / (double)RAND_MAX;
    // }



    s2 = clock();

    cpu1 = ((double)(s2 - s1)) / CLOCKS_PER_SEC;

    // Allocate memory for output arrays
    int* indices = (int*)malloc(num_points * k * sizeof(int)); // Output indices
    double* distances = (double*)malloc(num_points * k * sizeof(double)); // Output distances

    // Call the knnsearch function
    knnsearch_parallel(C, Q, num_points, dimensions, k, indices, distances);


    s3 = clock();

    // Example output (print the first point's k-nearest neighbors)
    printf("Nearest neighbors for the first query point:\n");
    for (int i = 0; i < k; ++i) {
        printf("Index: %d, Distance: %f\n", indices[i], distances[i]);
    }

    cpu2 = ((double)(s3 - s2)) / CLOCKS_PER_SEC;
    printf("%f\n%f\n", cpu1, cpu2);

    // Free allocated memory
    free(C);
    free(Q);
    free(indices);
    free(distances);

    return 0;
}
//to compile you will need: gcc -fopenmp -I/usr/include/openblas -L/usr/lib/openblas -o e1 e1.c -lopenblas -lm
