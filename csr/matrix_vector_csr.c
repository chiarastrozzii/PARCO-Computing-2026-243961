#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h> //in the sequential code, this is needed just for timing
#include <string.h>

#define N_RUNS 100

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

void matrix_vector_mul(const Sparse_CSR* sparse_csr, const double* vec, double* res, int parallel){
    if (parallel == 1){
        #pragma omp parallel for schedule(runtime)
        for (size_t i=0; i<sparse_csr->n_rows; ++i){
            res[i] = 0;
            size_t nz_start = sparse_csr->row_ptrs[i];
            size_t nz_end = sparse_csr->row_ptrs[i+1];
        
            for (size_t j = nz_start; j < nz_end; ++j) {
                size_t col = sparse_csr->col_indices[j];
                double val = sparse_csr->values[j];
                res[i] += val * vec[col];
            }
        }
    }else if(parallel == 2){
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < sparse_csr->n_rows; ++i) {
            double sum = 0.0;
            size_t nz_start = sparse_csr->row_ptrs[i];
            size_t nz_end   = sparse_csr->row_ptrs[i+1];

            #pragma omp simd reduction(+:sum)
            for (size_t j = nz_start; j < nz_end; ++j) {
                size_t col = sparse_csr->col_indices[j];
                sum += sparse_csr->values[j] * vec[col];
            }

            res[i] = sum;
        }

    }else{
        for (size_t i=0; i<sparse_csr->n_rows; ++i){
            res[i] = 0;
            size_t nz_start = sparse_csr->row_ptrs[i];
            size_t nz_end = sparse_csr->row_ptrs[i+1];
        
            for (size_t j = nz_start; j < nz_end; ++j) {
                size_t col = sparse_csr->col_indices[j];
                double val = sparse_csr->values[j];
                res[i] += val * vec[col];
            }
        }
    }
    

}



int main(int argc, char *argv[]){
    if (argc < 3){
        fprintf(stderr, "usage:\n");
        fprintf(stderr, "  %s SEQ <matrix_file>\n", argv[0]);
        fprintf(stderr, "  %s PAR <matrix_file> <num_threads> <schedule_type> <chunck_size>\n", argv[0]);
        fprintf(stderr, "schedule types: static | dynamic | guided\n");
        return 1;
    }

    int is_parallel = 0;
    int num_threads = 1;
    int chunck_size = 0;
    omp_sched_t schedule = omp_sched_static;


    if (strcmp(argv[1], "SEQ") == 0){
        is_parallel = 0;
    }else if(strcmp(argv[1], "PAR") == 0){
        if (argc < 6) {
            fprintf(stderr, "error: Missing arguments for PAR mode.\n");
            fprintf(stderr, "usage: %s PAR <matrix_file> <num_threads> <schedule_type>\n", argv[0]);
            return 1;
        }
        is_parallel = 1;
        num_threads = atoi(argv[3]); //gets from the terminal the number of threads
        chunck_size = atoi(argv[5]);

        if(strcmp(argv[4], "static") == 0){
            schedule = omp_sched_static;
        }else if(strcmp(argv[4], "dynamic") == 0){
            schedule = omp_sched_dynamic;
        }else if(strcmp(argv[4], "guided") == 0){
            schedule = omp_sched_guided;
        }else if(strcmp(argv[4], "auto") == 0){
            schedule = omp_sched_auto;
        }else{
            fprintf(stderr, "unknown schedule type: %s\n", argv[4]);
            return 1;
        }

        omp_set_num_threads(num_threads);
        omp_set_schedule(schedule, chunck_size);

    }else if(strcmp(argv[1], "SIMD") == 0){
        is_parallel = 2;
        schedule = omp_sched_static;
        num_threads = atoi(argv[3]); //gets from the terminal the number of threads
        chunck_size = atoi(argv[5]);

        omp_set_num_threads(num_threads);
        omp_set_schedule(schedule, chunck_size);

    }
    else{
        fprintf(stderr, "invalid mode: %s\n", argv[1]);
    }



    FILE *f = fopen(argv[2], "r");
    if (!f){
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

    int n_rows, n_columns, n_nz;
    sscanf(line, "%d %d %d", &n_rows, &n_columns, &n_nz);

    int *row = malloc(n_nz * sizeof(int)); //arrays where i can store the matrix, we store it in the memory dynamically
    int *col = malloc(n_nz * sizeof(int)); //malloc because we just need to allocate memory since we fill it immediately after
    double *val = malloc(n_nz * sizeof(double));
    
    for (int i = 0; i < n_nz; ++i) {
        int r, c;
        double v;
        fscanf(f, "%d %d %lf", &r, &c, &v);
        row[i] = r - 1;   // convert 1-based to 0-based
        col[i] = c - 1;   // convert 1-based to 0-based
        val[i] = v;
    }

    Sparse_CSR sparse_csr;
    create_sparse_csr(n_rows, n_columns, n_nz, row, col, val, &sparse_csr);
    //print_sparse_csr(&sparse_csr);

    double* vec = malloc(n_columns * sizeof(double));
    random_vector(vec, n_columns);

    double* res = malloc(n_rows * sizeof(double));

    printf("matrix considered: %s\n", argv[2]);
    printf("the mode used was: %s\n", argv[1]);
    if (is_parallel == 1) {
        printf("number of threads used: %d | type of schedule: %s | chunck size considered: %d\n", num_threads, argv[4], chunck_size);
    }else if(is_parallel == 2){
        printf("number of threads used: %d | type of schedule: static | chunck size considered: %d\n", num_threads, chunck_size);
    }


    //to benchmark
    double times[N_RUNS];
    for(int run=0; run<N_RUNS; ++run){ //the loop does the matrix multiplication multiple times (10) and saves how much time it took
        double start = omp_get_wtime();
        matrix_vector_mul(&sparse_csr, vec, res, is_parallel);
        double end = omp_get_wtime();
        times[run] = (end - start) * 1000.0; //milliseconds
    }

    double final_time = 0.0;
    
    for (size_t i = 0; i<N_RUNS; ++i){
        final_time+=times[i];
    }
    printf("average time took: %.3f ms\n", final_time/N_RUNS);
    printf("\n");

    free(row);
    free(col);
    free(val);
    free(vec);
    free(res);
    free_sparse_csr(&sparse_csr);

    return 0;

}