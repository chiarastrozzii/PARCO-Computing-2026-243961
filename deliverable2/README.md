# Deliverable 2 — Distributed SpMV (MPI 1D/2D) + Hybrid MPI/OpenMP

This deliverable implements and benchmarks **Sparse Matrix–Vector Multiplication (SpMV)** in **CSR format**, using:
- **MPI 1D cyclic row partitioning** with ghost/halo exchange
- **MPI 2D block partitioning** on a Cartesian process grid (p × q)
- Optional **OpenMP parallelization** for the local SpMV kernel (hybrid runs)

The executable supports both decomposition strategies and two execution modes (sequential local kernel vs OpenMP-parallel local kernel).

---

## Project Structure (Deliverable 2)

Inside the main repository, the root folder contains:
- `deliverable1/`
- `deliverable2/` ← this folder
- `matrix_mkt/` (shared matrices used across deliverables)

Within `deliverable2/`:

- `SpMV_MPI.c` — main program (MPI setup, distribution, timings, verification)
- `communications/` — communication utilities (scatter, gather, ghost exchange, 2D vector distribution)
- `multiplication/` — SpMV kernel (CSR) + OpenMP support
- `config/` — Matrix Market reader and helpers (CSR creation, utilities)

Data folders:
- `strong_scaling_mtx/` — real matrices (Matrix Market format) for strong scaling tests
- `random_matrices/` — synthetic matrices generated for weak scaling experiments

Scripts:
- `scripts/`
  - build scripts: one for **Linux** and one for **macOS**
  - benchmarking scripts:
    - strong scaling
    - weak scaling (MPI-only)
    - weak scaling (hybrid MPI+OpenMP)
  - Python generator for synthetic matrices

Cluster jobs:
- `pbs/` — PBS job files for the 3 benchmarking configurations

---

## Build Instructions (CMake)

### Requirements
- CMake ≥ 3.16
- MPI (C)
- OpenMP

### Build (MacOs)
```bash
cd scripts
./build_mac.sh 
```

### Build (Linux)
```bash
cd scripts
./build_linux.sh 
```

### To run a desired configuration
```bash
#from deliverable2 directory
cd build
mpirun -np <num_process> ./spmv <matrix> <1D/2D> <SEQ/PAR>
```

### To run scripts for scaling configurations
```bash
#from deliverable2 directory
cd build
../scripts/strong_scaling.sh #modify with wanted configuration (weak_scaling.sh or weak_hybrid.sh)
```


