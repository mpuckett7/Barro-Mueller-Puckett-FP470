#!/bin/bash
# Run to view full or partial results for MPI testing of the bifurcation eulerLagrange example

    echo "== Bifurcation, LaGrange with 1 processes, with contributions =="
    cat bifurcation3d1-1.txt
    echo

for n in 2 4 8 16 32; do
    echo "== Bifurcation, LaGrange with $n processes, with contributions =="
    cat bifurcation3d_changed-$n.txt
    echo
done