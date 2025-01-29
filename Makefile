# Makefile for both omp and cilk versions

# Compiler settings for each directory
OMP_CC = gcc

# Find OpenCilk clang by searching the opencilk directory pattern
CILK_CC = ~/opencilk/opencilk-2.1.0-x86_64-linux-gnu-ubuntu-22.04/bin/clang

# If CILK_CC is not set by the find, fall back to a default clang
CILK_CC ?= $(shell which clang)

# Compiler flags
OMP_FLAGS = -fopenmp -I /usr/include/openblas -L /usr/lib/openblas
CILK_FLAGS = -I /usr/include/openblas -L /usr/lib/openblas -fopencilk -O3

# Directories
OMP_DIR = omp_pthreads
CILK_DIR = cilk_pthreads

# Output filenames
OMP_OUTPUT = $(OMP_DIR)/omp_program
CILK_OUTPUT = $(CILK_DIR)/cilk_program

# Source files
OMP_SRC = $(OMP_DIR)/omp_pthreads.c
CILK_SRC = $(CILK_DIR)/cilk_pthreads.c

# Header files
OMP_HEADER = $(OMP_DIR)/knn_omp_pthreads.h
CILK_HEADER = $(CILK_DIR)/knn_cilk_pthreads.h

# Libraries
LIBS = -lopenblas -lm

# Targets
all: $(OMP_OUTPUT) $(CILK_OUTPUT)

# Rule to compile omp_program and place it in the omp_pthreads directory
$(OMP_OUTPUT): $(OMP_SRC) $(OMP_HEADER)
	$(OMP_CC) $(OMP_FLAGS) -o $(OMP_OUTPUT) $(OMP_SRC) $(LIBS)

# Rule to compile cilk_program and place it in the cilk_pthreads directory
$(CILK_OUTPUT): $(CILK_SRC) $(CILK_HEADER)
	$(CILK_CC) $(CILK_FLAGS) -o $(CILK_OUTPUT) $(CILK_SRC) $(LIBS)

# Clean up generated files
clean:
	rm -f $(OMP_OUTPUT) $(CILK_OUTPUT)

.PHONY: all clean
