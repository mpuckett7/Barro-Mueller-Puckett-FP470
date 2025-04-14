#!/bin/bash
# Run to view full or partial results for OpenMP testing of the EXAMPLE-NAME example


for n in 1 2  4 8 16 32 64; do
    echo "== Cylinder3d with $n threads, NO changes =="
    cat cylinder3d_unchanged-$n.txt
    echo
done

for n in 1 2 4 8 16 32 64; do
    echo "== Cylinder3d with $n threads, with contributions =="
    cat cylinder3d_changed-$n.txt
    echo
done