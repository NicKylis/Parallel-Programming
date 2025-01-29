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

## How to Use the Makefile

This project contains a `Makefile` to compile two different versions of a program: one using **OpenMP** and the other using **OpenCilk**. The Makefile automates the process of compiling both versions, each with the appropriate flags and compilers.

### Prerequisites

Before running the `Makefile`, make sure you have the following installed:

- **GCC**: For the OpenMP version.
- **OpenCilk**: For the OpenCilk version of the program.
  - You must have OpenCilk installed on your system, and you should be using the OpenCilk-enabled `clang` compiler (not the system default `clang`).
  
  **To install OpenCilk**, follow the instructions on the [OpenCilk GitHub page](https://github.com/OpenCilk).
  
  If you already have OpenCilk installed, make sure the `clang` compiler provided by OpenCilk is available in your system `PATH` or directly referenced in the `Makefile`.

- **OpenBLAS**: Make sure OpenBLAS is installed and the required headers and libraries are available:
  - **Headers**: `/usr/include/openblas`
  - **Libraries**: `/usr/lib/openblas`

### Quick Start

1. **Clone the repository** (if you haven't already):

   ```bash
   git clone https://github.com/NicKylis/Parallel-Programming.git
   cd Parallel-Programming
   ```

2. **Check that your `clang` points to the OpenCilk version**:
   
   Run the following to ensure the OpenCilk `clang` is being used:

   ```bash
   which clang
   clang --version
   ```

   If the output doesn't mention "OpenCilk," make sure to update your `PATH` to include the OpenCilk `clang`:

   ```bash
   export PATH=~/opencilk/opencilk-YOUR_VERSION_HERE/bin:$PATH
   ```

   To make this change permanent, add the `export PATH=...` line to your `~/.bashrc` (or `~/.zshrc` if you're using Zsh) and run:

   ```bash
   source ~/.bashrc
   ```

3. **Run the `Makefile`** to compile both versions of the program:

   ```bash
   make
   ```

   This will compile two executables:
   - **`omp_program`** in the `omp_pthreads` directory (compiled with OpenMP).
   - **`cilk_program`** in the `cilk_pthreads` directory (compiled with OpenCilk).

4. **Clean up compiled files**:
   
   To remove the compiled executables, run:

   ```bash
   make clean
   ```

   This will delete both `omp_program` and `cilk_program`.

### Troubleshooting

- If you encounter an error like `unknown argument: '-fopencilk'`, it means you're using the wrong cilk compiler.

## AUTHOR
Nikolaos Kylintireas