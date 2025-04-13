#!/bin/bash
# Run to view full or partial results for OpenMP testing of the bifurcation eulerLagrange example


for n in 1 2  4 8 16 32 64; do
    echo "== Bifurcation, LaGrange with $n threads, NO changes =="
    cat bifurcation3d_unchanged-$n.txt
    echo
done

for n in 1 2 4 8 16 32 64; do
    echo "== Bifurcation, LaGrange with $n threads, with contributions =="
    cat bifurcation3d_changed-$n.txt
    echo
done