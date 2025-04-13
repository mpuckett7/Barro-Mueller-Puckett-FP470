#!/bin/bash
# Run to view full or partial results for MPI testing of the bifurcation eulerLagrange example


for n in 1 2  4 8 16 32 64; do
    echo "== Bifurcation, LaGrange with $n processes, NO changes =="
    cat bifurcation3d_unchanged-$n.txt
    echo
done

for n in 1 2 4 8 16 32 64; do
    echo "== Bifurcation, LaGrange with $n processes, with contributions =="
    cat bifurcation3d_changed-$n.txt
    echo
done