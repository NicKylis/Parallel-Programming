#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <cblas.h>
#include <pthread.h>

#define min(a , b) ((a) < (b) ? (a) : (b))
#define NUM_THREADS 8 //Change this to the max cores of your system!!

// Structure to pass data to threads
typedef struct {
    double *distances;
    int *indices;
    int low;
    int count;
    int dir; // Sorting direction: 1 for ascending, 0 for descending, could be bool
} threadData;

// Helper function for quickselect
int partition(double* distances, int* indices, int left, int right, int pivot_idx);

// Sorting algorithm quickselect for knnsearch, as instructed
void quickselect(double* distances, int* indices, int left, int right, int k);

// A simple knnsearch algorithm divided by memory blocks (V0)
void knnsearch(double* C, double* Q, int num_points, int dimensions, int k, int* indices, double* distances);

// Similar to knnsearch, but using the openBLAS library for SOME of the arithmetic operations (seemed to be faster this way)
void knnsearchBLAS(double* C, double* Q, int num_points, int dimensions, int k, int* indices, double* distances);

// Parallel execution of the knnsearch, with a simple sorting algorithm
void knnsearch_parallel(double* C, double* Q, int num_points, int dimensions, int k, int* indices, double* distances);

// Helper functions for paralllel bitonic sort
void swap(double *a, double *b);
void swapIndices(int *a, int *b);
void compareAndSwap(double distances[], int indices[], int i, int j, int dir);
void bitonicMerge(double distances[], int indices[], int low, int count, int dir);
void *parallelBitonicSort(void *arg);
void bitonicSort(double distances[], int indices[], int n, int ascending);

// Similar to knnsearch_parallel but with parallel bitonic sorting for optimized performance (not working properly yet)
void knnsearch_parallelBITONIC(double* C, double* Q, int num_points, int dimensions, int k, int* indices, double* distances);

int partition(double* distances, int* indices, int left, int right, int pivot_idx) {
    double pivot_val = distances[pivot_idx];
    double temp = indices[pivot_idx];
    indices[pivot_idx] = indices[right];
    indices[right] = temp;
    
    double temp_dist = distances[pivot_idx];
    distances[pivot_idx] = distances[right];
    distances[right] = temp_dist;
    
    int store_idx = left;
    for (int i = left; i < right; i++) {
        if (distances[i] < pivot_val) {
            temp = indices[i];
            indices[i] = indices[store_idx];
            indices[store_idx] = temp;
            
            temp_dist = distances[i];
            distances[i] = distances[store_idx];
            distances[store_idx] = temp_dist;
            
            store_idx++;
        }
    }
    temp = indices[store_idx];
    indices[store_idx] = indices[right];
    indices[right] = temp;
    
    temp_dist = distances[store_idx];
    distances[store_idx] = distances[right];
    distances[right] = temp_dist;
    
    return store_idx;
}

void quickselect(double* distances, int* indices, int left, int right, int k) {
    if (left == right) return;
    
    int pivot_idx = left + rand() % (right - left + 1);
    pivot_idx = partition(distances, indices, left, right, pivot_idx);
    
    if (k == pivot_idx) return;
    else if (k < pivot_idx) quickselect(distances, indices, left, pivot_idx - 1, k);
    else quickselect(distances, indices, pivot_idx + 1, right, k);
}

//Here is my thought process, leaving this in case someone else wants to test the openBLAS library and how it could work faster

// void knnsearch(double* C, double* Q, int num_points, int dimensions, int k, int* indices, double* distances) {
//     int block_size = 100; // Define block size according to memory constraints
//     int num_blocks = (num_points + block_size - 1) / block_size;
    
//     #pragma omp parallel for
//     for (int i = 0; i < num_blocks; ++i) {
//         int start_idx = i * block_size;
//         int end_idx = min((i + 1) * block_size, num_points);
        
//         for (int j = start_idx; j < end_idx; ++j) {
//             double* dist_block = (double*)malloc(num_points * sizeof(double));

//             for (int l = 0; l < num_points; ++l) {
//                 double dist = 0.0;
                
//                 // Compute squared distance between point j in Q and point l in C
//                 for (int d = 0; d < dimensions; ++d) {
//                     double diff = Q[j * dimensions + d] - C[l * dimensions + d];
//                     dist += diff * diff;
//                 }
                
//                 dist_block[l] = sqrt(dist); // Element-wise sqrt
//             }




//             // double alpha = -2.0, beta = 0.0;
//             // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
//             //             num_points, 1, dimensions, alpha, C, dimensions, 
//             //             &Q[j * dimensions], dimensions, beta, dist_block, 1);
            
//             // // Calculate the square terms C^2 and Q^2 and add them to the distance block
//             // for (int l = 0; l < num_points; ++l) {
//             //     double c_sq = 0.0, q_sq = 0.0;
//             //     for (int d = 0; d < dimensions; ++d) {
//             //         c_sq += C[l * dimensions + d] * C[l * dimensions + d];
//             //         q_sq += Q[j * dimensions + d] * Q[j * dimensions + d];
//             //     }
//             //     dist_block[l] += c_sq + q_sq;
//             // }
            
//             // Apply element-wise square root to finalize distances
//             for (int l = 0; l < num_points; ++l) {
//                 dist_block[l] = sqrt(dist_block[l]);
//             }

//             // Find k smallest distances
//             int* idx_block = (int*)malloc(num_points * sizeof(int));
//             for (int m = 0; m < num_points; m++) idx_block[m] = m;
            
//             quickselect(dist_block, idx_block, 0, num_points - 1, k);
            
//             for (int m = 0; m < k; m++) {
//                 indices[j * k + m] = idx_block[m];
//                 distances[j * k + m] = dist_block[idx_block[m]];
//             }
            
//             free(dist_block);
//             free(idx_block);
//         }
//     }
// }

void knnsearch(double* C, double* Q, int num_points, int dimensions, int k, int* indices, double* distances) {
    int block_size = 100; // Define block size according to memory constraints
    int num_blocks = (num_points + block_size - 1) / block_size;

    #pragma omp parallel for
    for (int i = 0; i < num_blocks; ++i) {
        int start_idx = i * block_size;
        int end_idx = fmin((i + 1) * block_size, num_points);

        // Process each query point within this block
        for (int j = start_idx; j < end_idx; ++j) {
            // Allocate within the thread's scope to avoid shared memory issues
            double* dist_block = (double*)malloc(num_points * sizeof(double));
            int* idx_block = (int*)malloc(num_points * sizeof(int));

            // Calculate distances from query point Q[j] to each point in C
            for (int l = 0; l < num_points; ++l) {
                double dist = 0.0;

                // Compute squared distance between point j in Q and point l in C
                for (int d = 0; d < dimensions; ++d) {
                    double diff = Q[j * dimensions + d] - C[l * dimensions + d];
                    dist += diff * diff;
                }

                dist_block[l] = dist;  // Only square root after finding nearest neighbors
                idx_block[l] = l;      // Initialize index array
            }

            // Find k smallest distances using quickselect
            quickselect(dist_block, idx_block, 0, num_points - 1, k);

            // After quickselect, finalize k nearest distances by taking the square root
            for (int m = 0; m < k; ++m) {
                indices[j * k + m] = idx_block[m];
                distances[j * k + m] = sqrt(dist_block[idx_block[m]]);
            }

            // Free memory allocated within each thread
            free(dist_block);
            free(idx_block);
        }
    }
}

void knnsearchBLAS(double* C, double* Q, int num_points, int dimensions, int k, int* indices, double* distances) {
    int block_size = 100;
    int num_blocks = (num_points + block_size - 1) / block_size;

    // Precompute squared norms of each point in C
    double* C_squared_norms = (double*)malloc(num_points * sizeof(double));
    for (int i = 0; i < num_points; ++i) {
        double norm_sq = 0.0;
        for (int d = 0; d < dimensions; ++d) {
            norm_sq += C[i * dimensions + d] * C[i * dimensions + d];
        }
        C_squared_norms[i] = norm_sq;
    }

    #pragma omp parallel for
    for (int i = 0; i < num_blocks; ++i) {
        int start_idx = i * block_size;
        int end_idx = fmin((i + 1) * block_size, num_points);

        for (int j = start_idx; j < end_idx; ++j) {
            // Allocate distance and index storage for the current query point
            double* dist_block = (double*)malloc(num_points * sizeof(double));
            int* idx_block = (int*)malloc(num_points * sizeof(int));

            // Compute squared norm of the current query point Q[j]
            double Q_norm_sq = 0.0;
            for (int d = 0; d < dimensions; ++d) {
                Q_norm_sq += Q[j * dimensions + d] * Q[j * dimensions + d];
            }

            // Calculate -2 * Q[j] * C^T using BLAS
            double* cross_terms = (double*)malloc(num_points * sizeof(double));
            cblas_dgemv(CblasRowMajor, CblasNoTrans, num_points, dimensions,
                        -2.0, C, dimensions, &Q[j * dimensions], 1, 0.0, cross_terms, 1);

            // Compute the full squared distances for each point in C
            for (int l = 0; l < num_points; ++l) {
                dist_block[l] = Q_norm_sq + C_squared_norms[l] + cross_terms[l];
                idx_block[l] = l;  // Initialize index array
            }

            // Find k smallest distances using quickselect
            quickselect(dist_block, idx_block, 0, num_points - 1, k);

            // Store sorted k nearest neighbors with final sqrt distances
            for (int m = 0; m < k; m++) {
                indices[j * k + m] = idx_block[m];
                distances[j * k + m] = sqrt(dist_block[idx_block[m]]);
            }

            // Free memory allocated within each thread
            free(dist_block);
            free(idx_block);
            free(cross_terms);
        }
    }

    // Free memory allocated for squared norms of C
    free(C_squared_norms);
}


void knnsearch_parallel(double* C, double* Q, int num_points, int dimensions, int k, int* indices, double* distances) {
    // Split the dataset in half
    int half_points = num_points / 2;
    
    // Allocate memory for the indices and distances from the two halves
    int* indices_half1 = (int*)malloc(num_points * k * sizeof(int));
    int* indices_half2 = (int*)malloc(num_points * k * sizeof(int));
    double* distances_half1 = (double*)malloc(num_points * k * sizeof(double));
    double* distances_half2 = (double*)malloc(num_points * k * sizeof(double));

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            // Perform knnsearch on the first half
            knnsearchBLAS(C, Q, half_points, dimensions, k, indices_half1, distances_half1);
        }

        #pragma omp section
        {
            // Perform knnsearch on the second half
            knnsearchBLAS(C + half_points * dimensions, Q, half_points, dimensions, k, indices_half2, distances_half2);
        }
    }

    // Combine the results from the two halves
    for (int j = 0; j < num_points; ++j) {
        // Merge results for point j (from both halves)
        double* combined_distances = (double*)malloc(2 * k * sizeof(double));
        int* combined_indices = (int*)malloc(2 * k * sizeof(int));

        // Copy distances and indices from both halves
        for (int m = 0; m < k; ++m) {
            combined_distances[m] = distances_half1[j * k + m];
            combined_indices[m] = indices_half1[j * k + m];
            combined_distances[k + m] = distances_half2[j * k + m];
            combined_indices[k + m] = indices_half2[j * k + m];
        }

        // Sort the combined results by distance (using a simple sort algorithm)
        for (int m = 0; m < 2 * k - 1; ++m) {
            for (int n = m + 1; n < 2 * k; ++n) {
                if (combined_distances[m] > combined_distances[n]) {
                    // Swap distances
                    double temp_dist = combined_distances[m];
                    combined_distances[m] = combined_distances[n];
                    combined_distances[n] = temp_dist;

                    // Swap corresponding indices
                    int temp_idx = combined_indices[m];
                    combined_indices[m] = combined_indices[n];
                    combined_indices[n] = temp_idx;
                }
            }
        }

        // Store the final top-k indices and distances
        for (int m = 0; m < k; ++m) {
            indices[j * k + m] = combined_indices[m];
            distances[j * k + m] = combined_distances[m];
        }

        // Free allocated memory
        free(combined_distances);
        free(combined_indices);
    }

    // Free the memory for half results
    free(indices_half1);
    free(indices_half2);
    free(distances_half1);
    free(distances_half2);
}

// Swap function for both distances and indices
void swap(double *a, double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

void swapIndices(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Compare and Swap function
void compareAndSwap(double distances[], int indices[], int i, int j, int dir) {
    if ((dir == 1 && distances[i] > distances[j]) || (dir == 0 && distances[i] < distances[j])) {
        // Swap distances
        swap(&distances[i], &distances[j]);
        // Swap corresponding indices
        swapIndices(&indices[i], &indices[j]);
    }
}

// Merge function for bitonic sort
void bitonicMerge(double distances[], int indices[], int low, int count, int dir) {
    if (count > 1) {
        int k = count / 2;
        for (int i = low; i < low + k; i++) {
            compareAndSwap(distances, indices, i, i + k, dir);
        }
        bitonicMerge(distances, indices, low, k, dir);
        bitonicMerge(distances, indices, low + k, k, dir);
    }
}

// Thread function for recursive parallel bitonic sort
void *parallelBitonicSort(void *arg) {
    threadData *data = (threadData *)arg;
    double *distances = data->distances;
    int *indices = data->indices;
    int low = data->low;
    int count = data->count;
    int dir = data->dir;

    if (count > 1) {
        int k = count / 2;

        // Prepare thread data for ascending and descending halves
        threadData left = {distances, indices, low, k, 1};
        threadData right = {distances, indices, low + k, k, 0};

        pthread_t thread1, thread2;

        // Create threads for the two recursive sorts
        pthread_create(&thread1, NULL, parallelBitonicSort, &left);
        pthread_create(&thread2, NULL, parallelBitonicSort, &right);

        // Wait for both threads to complete
        pthread_join(thread1, NULL);
        pthread_join(thread2, NULL);

        // Merge the sorted sequences
        bitonicMerge(distances, indices, low, count, dir);
    }

    return NULL;
}

// Entry point for bitonic sort
void bitonicSort(double distances[], int indices[], int n, int ascending) {
    threadData data = {distances, indices, 0, n, ascending};
    parallelBitonicSort(&data);
}

// Main k-NN search function
void knnsearch_parallelBITONIC(double* C, double* Q, int num_points, int dimensions, int k, int* indices, double* distances) {
    // Split the dataset in half
    int half_points = num_points / 2;

    // Allocate memory for the indices and distances from the two halves
    int* indices_half1 = (int*)malloc(num_points * k * sizeof(int));
    int* indices_half2 = (int*)malloc(num_points * k * sizeof(int));
    double* distances_half1 = (double*)malloc(num_points * k * sizeof(double));
    double* distances_half2 = (double*)malloc(num_points * k * sizeof(double));

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            // Perform knnsearch on the first half
            knnsearchBLAS(C, Q, half_points, dimensions, k, indices_half1, distances_half1);
        }

        #pragma omp section
        {
            // Perform knnsearch on the second half
            knnsearchBLAS(C + half_points * dimensions, Q, half_points, dimensions, k, indices_half2, distances_half2);
        }
    }

    // Combine the results from the two halves
    for (int j = 0; j < num_points; ++j) {
        // Merge results for point j (from both halves)
        double* combined_distances = (double*)malloc(2 * k * sizeof(double));
        int* combined_indices = (int*)malloc(2 * k * sizeof(int));

        // Copy distances and indices from both halves
        for (int m = 0; m < k; ++m) {
            combined_distances[m] = distances_half1[j * k + m];
            combined_indices[m] = indices_half1[j * k + m];
            combined_distances[k + m] = distances_half2[j * k + m];
            combined_indices[k + m] = indices_half2[j * k + m];
        }

        // Sort the combined results by distance using parallel bitonic sort
        bitonicSort(combined_distances, combined_indices, 2 * k, 1); // Sort in ascending order of distance

        // Store the final top-k indices and distances
        for (int m = 0; m < k; ++m) {
            indices[j * k + m] = combined_indices[m];
            distances[j * k + m] = combined_distances[m];
        }

        // Free allocated memory for combined results
        free(combined_distances);
        free(combined_indices);
    }

    // Free the memory for half results
    free(indices_half1);
    free(indices_half2);
    free(distances_half1);
    free(distances_half2);
}

#endif