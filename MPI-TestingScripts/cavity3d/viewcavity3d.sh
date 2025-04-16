#!/bin/bash
# Run to view full or partial results for MPI testing of the cavity3d example

for n in 1 2 4 8 16 32 64; do
    echo "== cavity3d with $n processes =="
    cat cavity3d_changed-$n.txt
    echo
done