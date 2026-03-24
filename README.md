# PARCO Computing — Project Deliverables

This repository contains the implementations and experiments for the PARCO Computing project on **Sparse Matrix–Vector Multiplication (SpMV)** using OpenMP and MPI.

## Repository Structure

```text
deliverables/
├── deliverable1/        # Parallel SpMV in CSR format parallelized using OpenMP
├── deliverable2/        # Distributed SpMV with 1D/2D MPI decomposition and hybrid MPI+OpenMP
├── matrix_mkt/          # Shared Matrix Market matrices used across deliverables for small comparison
└── README.md            # Project navigation
```

## Deliverables

### 📁 Deliverable 1
Study on the performance of the SpMV in CSR format with different OpenMP scheduling policies employed.

→ See `deliverable1/README.md` for build instructions and usage.

→ Technical Report: [Analysis of OpenMP Scheduling (PDF)](reports/Strozzi-D1-243961.pdf)

### 📁 Deliverable 2
Contains the distributed and hybrid parallel implementations:
- MPI 1D cyclic row partitioning
- MPI 2D block decomposition
- Hybrid MPI + OpenMP execution
- Strong and weak scaling experiments

→ See `deliverable2/README.md` for detailed documentation.

→ Technical Report: [Analysis of MPI Scaling (PDF)](reports/Strozzi-243961-D2.pdf)

---

Each deliverable is self-contained with its own build system, scripts, and experiment setup.
