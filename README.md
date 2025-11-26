# Parallel Sparse Matrix-Vector Multiplication with OpenMP

## 📌 Overview
This project evaluates the performance of **Sparse Matrix-Vector Multiplication (SpMV)** implemented using the **Compressed Sparse Row (CSR)** format and parallelized with **OpenMP**.

The project includes:

- **Sequential** SpMV implementation  
- **Parallel** SpMV implementation with OpenMP scheduling policies: `static`, `dynamic`, `guided`, `auto`, `runtime`  
- **SIMD-optimized** inner loop version  
- **Block-based CSB (Compressed Sparse Blocks)** variant  
- Automated scripts to run experiments, collect timings, and generate CSV output  

Performance was tested on multiple matrices with different sparsity patterns and thread counts up to 64.
