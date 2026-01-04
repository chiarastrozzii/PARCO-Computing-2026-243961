#ifndef CONFIG_H
#define CONFIG_H

#include <stddef.h>

typedef struct Sparse_CSR{
    size_t n_rows;
    size_t n_cols;
    size_t n_nz;

    size_t* row_ptrs;
    size_t* col_indices;
    double* values;
} Sparse_CSR;

void create_sparse_csr(
    size_t n_rows,
    size_t n_cols,
    size_t n_nz,

    const int *row,
    const int* col,
    const double* val,
    Sparse_CSR* output_csr
);
void print_sparse_csr(Sparse_CSR* sparse_csr);
void free_sparse_csr(Sparse_CSR* sparse_csr);
void read_matrix_market_file(
    const char* filename,
    int* n_rows,
    int* n_cols,
    int* n_nz,
    int** row_indices,
    int** col_indices,
    double** values
);
void random_vector(double* vec, size_t size);

#endif // CONFIG_H