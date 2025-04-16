#!/bin/bash
# Run to view full or partial results for OpenMP testing of the cavity3d example


for n in 1 2  4 8 16 32 64; do
    echo "== Cavity3d with $n threads, NO changes =="
    cat cavity3d_unchanged-$n.txt
    echo
done

for n in 1 2 4 8 16 32 64; do
    echo "== Cavity3d with $n threads, with contributions =="
    cat cavity3d_changed-$n.txt
    echo
done