# Parallel-Programming
The first excercise of the course of Parallel Programming. The objective is to create a knnsearch function in c and optimize it using parallel programming.
The project is divided in 4 directories:
1. cilk_pthreads
    This directory contains the implementation of both simple and parallel knnsearch using the OpenCilk and Pthreads libraries. Both a source file (.c) and
    a header file (.h) are included and needed. The source file contains the main function while the header one has all the other functions.
    To compile it in a terminal, the following command should be used, with the necessary adjustments for your system:
    /yourPath/clang -I/usr/include/openblas -L/usr/lib/openblas -o cilk_pthreads cilk_pthreads.c -lopenblas -lm -fopencilk -O3
2. omp_pthreads
    Similarly, this directory contains the implementation of both simple and parallel knnsearch using the OpenMP and Pthreads libraries. To compile it in a terminal,
    the following command should be used, with the necessary adjustments for your system:
    gcc -fopenmp -I/usr/include/openblas -L/usr/lib/openblas -o omp_pthreads omp_pthreads.c -lopenblas -lm
3. test_files
    This directory was created mostly for testing and having a natural flaw when coding. It reveals the whole process of the project
4. txtDataset
    This directory contains all the .txt files used, although for your convenience they have already been placed inside the directories.