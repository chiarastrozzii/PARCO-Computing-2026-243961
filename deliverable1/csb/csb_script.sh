#!/bin/bash

FILE=$1
OUTPUT="report_CSB.txt"

MATRICES=("matrix_mkt/1138_bus.mtx" "matrix_mkt/494_bus.mtx" "matrix_mkt/bcspwr06.mtx" "matrix_mkt/bcspwr08.mtx" "matrix_mkt/bcspwr10.mtx")
SCHEDULES=("static" "dynamic" "guided" "auto")
THREADS=(64 32 16 4)
CHUNKS=(0 1 10 50 100)

echo "Compiling $FILE with OpenMP..."
gcc-15 -fopenmp "$FILE" -o sparse_seq
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi
echo "Compilation successful."
echo "" > "$OUTPUT"

echo "Running sequential CSB..."
for MATRIX in "${MATRICES[@]}"; do
    ./sparse_seq SEQ_CSB "$MATRIX" >> "$OUTPUT"
done
echo "" >> "$OUTPUT"

echo "Running parallel CSB..."
for MATRIX in "${MATRICES[@]}"; do
    for SCHED in "${SCHEDULES[@]}"; do
        for NTHREAD in "${THREADS[@]}"; do
            for CHUNK in "${CHUNKS[@]}"; do
                export OMP_NUM_THREADS=$NTHREAD
                ./sparse_seq PAR "$MATRIX" $NTHREAD $SCHED $CHUNK >> "$OUTPUT"
            done
        done
    done
done
echo "" >> "$OUTPUT"

echo "All runs finished. Results saved to $OUTPUT."