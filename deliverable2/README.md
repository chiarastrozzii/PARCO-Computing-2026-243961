# Deliverable 2 — Distributed SpMV (MPI 1D/2D) + Hybrid MPI/OpenMP

This deliverable implements and benchmarks **Sparse Matrix–Vector Multiplication (SpMV)** in **CSR format**, using:
- **MPI 1D cyclic row partitioning** with ghost exchange
- **MPI 2D block partitioning** on a Cartesian process grid (p × q)
- Optional **OpenMP parallelization** for the local SpMV kernel (hybrid runs)

The executable supports both decomposition strategies and two execution modes (sequential local kernel vs OpenMP-parallel local kernel).
All experiments were perfomed on the institutional cluster provided with 142 CPU nodes (7674 cores), 10 GPU nodes and 65 TB of RAM.

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
mpirun -np <P> ./spmv <matrix_file> <1D/2D> <SEQ/PAR>
```
Arguments:
  - `<matrix_file>`: Matrix Market .mtx file
  - `<1D/2D>`: 
    - `1D`: cyclic row partitioning + ghost exchangew
    - `2D`: 2D block partitioning on Cartesian Grid
  - `<SEQ/PAR>`:
    - `SEQ` = local SpMV computed sequentially
    - `PAR` = OpenMP enabled
   
**Example PAR**
```bash
cd build
export OMP_NUM_THREADS=8
mpirun -np 16 ./spmv random_matrices/A_weak_1024.mtx 2D PAR
```

### To run scripts for scaling configurations
```bash
#from deliverable2 directory
./scripts/strong_scaling.sh #modify with wanted configuration (weak_scaling.sh or weak_hybrid.sh)
```
**Benchmarking scripts:**
- Strong scaling: fixed problem size, increasing MPI ranks
- Weak scaling (MPI) script: matrix size grows with MPI ranks
- Weak scaling (hybrid) script: matrix size grows with total threads (MPI ranks × OMP threads)

## Input Matrices

### Strong-scaling matrices (not included in the repository)
The following Matrix Market matrices were used for the strong scaling experiments:
- `af_2_k101.mtx`
- `af_shell10.mtx`
- `bone010.mtx`
- `dielFilterV2real.mtx`

These files are **not committed to the repository due to size constraints**.  
To reproduce the experiments, download them from the **SuiteSparse Matrix Collection** and place them in:
`deliverable2/strong_scaling_mtx/`

### Weak-scaling matrices
Weak-scaling experiments use **synthetic sparse matrices** generated with the provided Python script:
```bash
#from deliverable2 directory
cd scripts
python3 generator_random_mtx.py
```
The script outputs matrices in: `deliverable2/random_matrices/`

## PBS job submission
PBS jobs are provided in `/pbs` for:
- strong scaling
- weak scaling (MPY only)
- weak scaling (hybrid MPI + OpenMP)

Typical usage: 
```bash
qsub pbs/<job_file>.pbs
```

## Notes
When running the experiments on cluster some modifications are needed:
**In benchmarking scripts**
- remove building step as it needs to be added to `.pbs` file.
  ```bash
  #REMOVE:
  echo "Building project (cd .. && ./scripts/build_mac.sh) ..."
  pushd "$PROJECT_ROOT" >/dev/null
  ./scripts/build_mac.sh #REMOVE IN CLUSTER since build moved inside pbs file
  popd >/dev/null
  echo "Build done."
  echo
  ```
- complete number of processes up to 128 in both strong scaling and weak scaling
  - strong scaling -> `PROCS_LIST_DEFAULT=("1" "2" "4" "8" "16" "32" "64" "128")`
  - weak scaling (pure MPI) ->
    ```bash
    CONFIG_N=(1000 2000 4000 8000 16000 32000 64000 128000) #modify in cluster till$
    CONFIG_NNZ=(100000 200000 400000 800000 1600000 3200000 6400000 12800000)
    CONFIG_NP=(1 2 4 8 16 32 64 128)
    ```
  - weak scaling (OpenMP + MPI) -> `TOTAL_THREADS_LIST=("1" "2" "4" "8" "16" "32" "64" "128")`
    
**In `.pbs` files**:

Add this line before running the chosen configuration:
`./scripts/build_linux.sh`









