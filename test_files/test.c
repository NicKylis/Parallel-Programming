#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <cblas.h>
#include "functions.h"
#include <time.h>
#include <hdf5.h>
#include <pthread.h>

#define NUM_POINTS 10000 // Example number of points
#define DIMENSIONS 128   // Example dimensionality of each point
#define K 10             // Number of nearest neighbors

// int main() {
//     srand(100);
//     clock_t s1, s2, s3;
//     double cpu1, cpu2;
//    int num_points = NUM_POINTS;
//     int dimensions = DIMENSIONS;
//     int k = K;

//     // Allocate memory for C and Q matrices
//     double* C = (double*)malloc(num_points * dimensions * sizeof(double));
//     double* Q = (double*)malloc(num_points * dimensions * sizeof(double));

//     // Open the HDF5 file
//     const char *filename = "testfile.h5";
//     hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
//     if (file < 0) {
//         fprintf(stderr, "Error opening file %s\n", filename);
//         exit(0);
//         return EXIT_FAILURE;
//     }

//     // Open the datasets for C and Q matrices
//     hid_t dataset_C = H5Dopen(file, "/C", H5P_DEFAULT);  // Dataset C in the HDF5 file

//     if (dataset_C < 0) {
//         fprintf(stderr, "Error opening datasets in the file\n");
//         H5Fclose(file);
//         return EXIT_FAILURE;
//     }

//     // Read the C dataset into memory
//     herr_t status;
//     status = H5Dread(dataset_C, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, C);
//     if (status < 0) {
//         fprintf(stderr, "Error reading dataset C\n");
//         H5Dclose(dataset_C);
//         H5Fclose(file);
//         return EXIT_FAILURE;
//     }

//     H5Dclose(dataset_C);

//     // Close the HDF5 file
//     H5Fclose(file);

//     s1 = clock();

//     // Example initialization (of random data)
//     // for (int i = 0; i < num_points * dimensions; ++i) {
//     //     C[i] = rand() / (double)RAND_MAX;
//     //     Q[i] = rand() / (double)RAND_MAX;
//     // }

//     s2 = clock();

//     cpu1 = ((double)(s2 - s1)) / CLOCKS_PER_SEC;

//     // Allocate memory for output arrays
//     int* indices = (int*)malloc(num_points * k * sizeof(int)); // Output indices
//     double* distances = (double*)malloc(num_points * k * sizeof(double)); // Output distances

//     // Call the knnsearch function
//     knnsearch_parallelBITONIC(C, C, num_points, dimensions, k, indices, distances);


//     s3 = clock();

//     // Example output (print the first point's k-nearest neighbors)
//     printf("Nearest neighbors for the first query point:\n");
//     for (int i = 0; i < k; ++i) {
//         printf("Index: %d, Distance: %f\n", indices[i], distances[i]);
//     }

//     cpu2 = ((double)(s3 - s2)) / CLOCKS_PER_SEC;
//     printf("%f\n%f\n", cpu1, cpu2);

//     // Free allocated memory
//     free(C);
//     // free(Q);
//     free(indices);
//     free(distances);

//     return 0;
// }
//to compile you will need: gcc -fopenmp -I/usr/include/openblas -L/usr/lib/openblas -o test test.c -lopenblas -lm
//TODO optimize the parallel algorithm
//gcc -fopenmp -I/usr/include/openblas -L/usr/lib/openblas -o test test.c -lopenblas -lm -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5

int main() {
    srand(100);
    clock_t s1, s2, s3;
    double cpu1, cpu2;
    // Example parameters
    int num_points = 10000;      // Number of points in C and Q
    int dimensions = 32;         // Dimensionality of each point
    int k = 10;                  // Number of nearest neighbors to find

    // Allocate memory for C and Q
    double* C = (double*)malloc(num_points * dimensions * sizeof(double));
    double* Q = (double*)malloc(num_points * dimensions * sizeof(double));

    FILE *file = fopen("arrays.txt", "r"); //Write the filename.txt here
    
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

