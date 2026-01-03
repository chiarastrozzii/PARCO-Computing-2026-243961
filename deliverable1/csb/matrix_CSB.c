#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h> //in the sequential code, this is needed just for timing
#include <string.h>
#include <math.h>

typedef struct Sparse_CSR{
    size_t n_rows;
    size_t n_cols;
    size_t n_nz;

    size_t* row_ptrs;
    size_t* col_indices;
    double* values;
} Sparse_CSR;

typedef struct Sparse_CSB{
    int n_rows;
    int n_cols;
    int nnz;

    int beta;           
    int n_block_rows;   
    int n_block_cols;   
    int n_blocks;      

    int *blk_ptr;   
    int *row_ind;  
    int *col_indices;  
    double *val;
}Sparse_CSB;



void create_sparse_csr(
    size_t n_rows,
    size_t n_cols,
    size_t n_nz,

    const int *row,
    const int* col,
    const double* val,
    Sparse_CSR* output_csr
){
    output_csr->n_rows = n_rows;
    output_csr->n_cols = n_cols;
    output_csr->n_nz = n_nz;

    output_csr->row_ptrs = calloc(n_rows + 1, sizeof(size_t)); //calloc-> initializes memory and also fill it with zero, useful for initializing the pointer
    output_csr->col_indices = calloc(n_nz, sizeof(size_t));
    output_csr->values = calloc(n_nz, sizeof(double));

    //count non-zero per rows
    for (size_t i = 0; i < n_nz; ++i) {
        output_csr->row_ptrs[row[i] + 1]++;
    }

    for (size_t i = 0; i<n_rows; ++i){
        output_csr->row_ptrs[i + 1] += output_csr->row_ptrs[i];
    }

    //we fill the indeces of the columns and the values
    size_t *row_offset = calloc(n_rows, sizeof(size_t)); // track current position in each row
    for (size_t i = 0; i < n_nz; ++i) {
        int r = row[i];
        size_t dest = output_csr->row_ptrs[r] + row_offset[r];
        output_csr->col_indices[dest] = col[i];
        output_csr->values[dest] = val[i];
        row_offset[r]++;
    }

    free(row_offset);

}

void print_sparse_csr(Sparse_CSR* sparse_csr){
    printf("\n");
    printf("row\tcol\tval\n");
    printf("---\n");
    for (size_t i=0; i<sparse_csr->n_rows; ++i){
        size_t nz_start = sparse_csr->row_ptrs[i];
        size_t nz_end = sparse_csr->row_ptrs[i+1];
        
        for (size_t j = nz_start; j < nz_end; ++j) {
            size_t col = sparse_csr->col_indices[j];
            double val = sparse_csr->values[j];
            printf("%zu\t%zu\t%.2f\n", i, col, val);
        }
    }
}

void free_sparse_csr(Sparse_CSR* sparse_csr){
    free(sparse_csr->row_ptrs);
    free(sparse_csr-> col_indices);
    free(sparse_csr-> values);
}

void random_vector(double* vec, size_t size){
   static bool seeded = false;
    if (!seeded) { //sees run one, so that multiple calls for benchmar not reset the seed
        srand(time(NULL));
        seeded = true;
    }

    for(size_t i = 0; i < size; ++i){
        vec[i] = rand() % 10; 
    }
}

//converting the CSR matrix into CSB
static int iceil(int a, int b) {
    return (a + b - 1) / b;
}


static int compute_beta(int n_rows) {
    if (n_rows <= 0) return 1;
    int b = (int)floor(sqrt((double)n_rows));
    if (b < 1) b = 1;
    return b;
}


int csr_to_csb(const Sparse_CSR *csr, Sparse_CSB *csb)
{
    if (!csr || !csb) return -1;

    int n_rows = csr->n_rows;
    int n_cols = csr->n_cols;
    int nnz    = csr->n_nz;

    int beta = compute_beta(n_rows); 

    int n_block_rows = iceil(n_rows, beta);
    int n_block_cols = iceil(n_cols, beta);
    int n_blocks     = n_block_rows * n_block_cols;


    csb->n_rows       = n_rows;
    csb->n_cols       = n_cols;
    csb->nnz          = nnz;
    csb->beta         = beta;
    csb->n_block_rows = n_block_rows;
    csb->n_block_cols = n_block_cols;
    csb->n_blocks     = n_blocks;

    csb->blk_ptr = (int*)calloc(n_blocks + 1, sizeof(int));
    csb->row_ind = (int*)malloc(nnz * sizeof(int));
    csb->col_indices = (int*)malloc(nnz * sizeof(int));
    csb->val     = (double*)malloc(nnz * sizeof(double));

    if (!csb->blk_ptr || !csb->row_ind || !csb->col_indices || !csb->val) {
        free(csb->blk_ptr);
        free(csb->row_ind);
        free(csb->col_indices);
        free(csb->val);
        return -1;
    }


    int *block_nnz = (int*)calloc(n_blocks, sizeof(int));
    if (!block_nnz) {
        free(csb->blk_ptr);
        free(csb->row_ind);
        free(csb->col_indices);
        free(csb->val);
        return -1;
    }

    for (int i = 0; i < n_rows; ++i) {
        int row_start = csr->row_ptrs[i];
        int row_end   = csr->row_ptrs[i + 1];

        int block_row = i / beta;

        for (int k = row_start; k < row_end; ++k) {
            int j = csr->col_indices[k];
            int block_col = j / beta;
            int b = block_row * n_block_cols + block_col;
            block_nnz[b]++;
        }
    }


    csb->blk_ptr[0] = 0;
    for (int b = 0; b < n_blocks; ++b) {
        csb->blk_ptr[b + 1] = csb->blk_ptr[b] + block_nnz[b];
    }

    for (int b = 0; b < n_blocks; ++b) {
        block_nnz[b] = csb->blk_ptr[b];
    }

    for (int i = 0; i < n_rows; ++i) {
        int row_start = csr->row_ptrs[i];
        int row_end   = csr->row_ptrs[i + 1];

        int block_row = i / beta;
        int local_row = i % beta;

        for (int k = row_start; k < row_end; ++k) {
            int j = csr->col_indices[k];
            double v = csr->values[k];

            int block_col = j / beta;
            int local_col = j % beta;
            int b = block_row * n_block_cols + block_col;

            int pos = block_nnz[b]++; 

            csb->row_ind[pos] = local_row;
            csb->col_indices[pos] = local_col;
            csb->val[pos]     = v;
        }
    }

    free(block_nnz);
    return 0;
}

void free_csb(Sparse_CSB *csb)
{
    if (!csb) return;
    free(csb->blk_ptr);
    free(csb->row_ind);
    free(csb->col_indices);
    free(csb->val);

    csb->blk_ptr = NULL;
    csb->row_ind = NULL;
    csb->col_indices = NULL;
    csb->val     = NULL;
}

void csb_spmv_sequential(const Sparse_CSB *A, const double *x, double *y)
{
    int n_rows = A->n_rows;
    int n_cols = A->n_cols;
    int beta   = A->beta;
    int n_block_rows = A->n_block_rows;
    int n_block_cols = A->n_block_cols;


    for (int i = 0; i < n_rows; i++)
        y[i] = 0.0;

    for (int br = 0; br < n_block_rows; br++) {

        int row_start_global = br * beta;
        int row_end_global = row_start_global + beta;
        if (row_end_global > n_rows) row_end_global = n_rows;

        for (int bc = 0; bc < n_block_cols; bc++) {

            int block_id = br * n_block_cols + bc;
            int start = A->blk_ptr[block_id];
            int end   = A->blk_ptr[block_id + 1];

            int col_start_global = bc * beta;
            int col_end_global   = col_start_global + beta;
            if (col_end_global > n_cols) col_end_global = n_cols;

            for (int k = start; k < end; k++) {

                int local_r = A->row_ind[k];
                int local_c = A->col_indices[k];

                int gr = row_start_global + local_r;
                int gc = col_start_global + local_c;

                if (gr < n_rows && gc < n_cols)
                    y[gr] += A->val[k] * x[gc];
            }
        }
    }
}


void csb_spmv(const Sparse_CSB *A, const double *x, double *y)
{
    int n_rows       = A->n_rows;
    int n_cols       = A->n_cols;
    int beta         = A->beta;
    int n_block_rows = A->n_block_rows;
    int n_block_cols = A->n_block_cols;

    #pragma omp parallel for
    for (int i = 0; i < n_rows; ++i) {
        y[i] = 0.0;
    }

    #pragma omp parallel for schedule(dynamic)
    for (int br = 0; br < n_block_rows; ++br) {

        int row_start_global = br * beta;
        int row_end_global   = row_start_global + beta;
        if (row_end_global > n_rows) row_end_global = n_rows;

        for (int bc = 0; bc < n_block_cols; ++bc) {
            int block_id = br * n_block_cols + bc;
            int start = A->blk_ptr[block_id];
            int end   = A->blk_ptr[block_id + 1];

            int col_start_global = bc * beta;
            int col_end_global   = col_start_global + beta;
            if (col_end_global > n_cols) col_end_global = n_cols;

            for (int k = start; k < end; ++k) {
                int local_r = A->row_ind[k];
                int local_c = A->col_indices[k];

                int gr = row_start_global + local_r;
                int gc = col_start_global + local_c;

                if (gr < n_rows && gc < n_cols) {
                    y[gr] += A->val[k] * x[gc];
                }
            }
        }
    }
}




#define N_RUNS 10

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "usage:\n");
        fprintf(stderr, "  %s SEQ <matrix_file>\n", argv[0]);
        fprintf(stderr, "  %s PAR <matrix_file> <num_threads> <schedule_type> <chunk_size>\n", argv[0]);
        fprintf(stderr, "schedule types: static | dynamic | guided | auto\n");
        return 1;
    }

    bool is_parallel = false;
    bool seq_csb = false;
    int num_threads = 1;
    int chunk_size = 0;
    omp_sched_t schedule = omp_sched_static;

    // --- Parse mode and parameters ---
   if(strcmp(argv[1], "SEQ_CSB") == 0){
        is_parallel = false;
        seq_csb = true;
    } else if (strcmp(argv[1], "PAR") == 0) {
        if (argc < 6) {
            fprintf(stderr, "error: Missing arguments for PAR mode.\n");
            fprintf(stderr, "usage: %s PAR <matrix_file> <num_threads> <schedule_type> <chunk_size>\n", argv[0]);
            return 1;
        }
        is_parallel = true;
        num_threads = atoi(argv[3]);
        chunk_size = atoi(argv[5]);

        if (strcmp(argv[4], "static") == 0)
            schedule = omp_sched_static;
        else if (strcmp(argv[4], "dynamic") == 0)
            schedule = omp_sched_dynamic;
        else if (strcmp(argv[4], "guided") == 0)
            schedule = omp_sched_guided;
        else if (strcmp(argv[4], "auto") == 0)
            schedule = omp_sched_auto;
        else {
            fprintf(stderr, "unknown schedule type: %s\n", argv[4]);
            return 1;
        }

        omp_set_num_threads(num_threads);
        omp_set_schedule(schedule, chunk_size);
    } else {
        fprintf(stderr, "invalid mode: %s\n", argv[1]);
        return 1;
    }

    FILE *f = fopen(argv[2], "r");
    if (!f) {
        perror("error opening matrix file");
        return 1;
    }

    char line[256];
    do {
        if (!fgets(line, sizeof(line), f)) {
            fprintf(stderr, "Unexpected end of file\n");
            return 1;
        }
    } while (line[0] == '%');

    int n_rows, n_cols, n_nz;
    sscanf(line, "%d %d %d", &n_rows, &n_cols, &n_nz);

    int *row = malloc(n_nz * sizeof(int));
    int *col = malloc(n_nz * sizeof(int));
    double *val = malloc(n_nz * sizeof(double));

    for (int i = 0; i < n_nz; ++i) {
        int r, c;
        double v;
        fscanf(f, "%d %d %lf", &r, &c, &v);
        row[i] = r - 1;
        col[i] = c - 1;
        val[i] = v;
    }
    fclose(f);

    Sparse_CSR csr;
    create_sparse_csr(n_rows, n_cols, n_nz, row, col, val, &csr);

    Sparse_CSB csb;
    csr_to_csb(&csr, &csb);

    double *vec = malloc(n_cols * sizeof(double));
    double *res = malloc(n_rows * sizeof(double));
    random_vector(vec, n_cols);

    printf("matrix considered: %s\n", argv[2]);
    if (seq_csb){
        printf("the mode used was: sequential (CSB)\n");
    }else if (!is_parallel){
        printf("the mode used was: sequential (CSR)\n");
    }else{
        printf("the mode used was: parallel (CSB)\n");
    }
    

    

    if (is_parallel) {
        printf("number of threads used: %d | schedule: %s | chunk: %d\n",
               num_threads, argv[4], chunk_size);
    }

    double times[N_RUNS];
    for (int run = 0; run < N_RUNS; ++run) {
        double start = omp_get_wtime();

        if (is_parallel) {
            csb_spmv(&csb, vec, res);
        } else if(seq_csb){
            csb_spmv_sequential(&csb, vec, res);
        }
        double end = omp_get_wtime();
        times[run] = (end - start) * 1000.0; // ms
    }

    double avg = 0.0;
    for (int i = 0; i < N_RUNS; ++i) avg += times[i];
    avg /= N_RUNS;

    printf("average time took: %.3f ms\n\n", avg);

    free(row);
    free(col);
    free(val);
    free(vec);
    free(res);
    free_sparse_csr(&csr);
    free_csb(&csb);

    return 0;
}