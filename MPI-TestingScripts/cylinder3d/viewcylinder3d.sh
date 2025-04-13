#!/bin/bash
# Run to view full or partial results for MPI testing of the cylinder3d example


for n in 1 2  4 8 16 32 64; do
    echo "== cylinder3d with $n processes, NO changes =="
    cat cylinder3d_unchanged-$n.txt
    echo
done

for n in 1 2 4 8 16 32 64; do
    echo "== cylinder3d with $n processes, with contributions =="
    cat cylinder3d_changed-$n.txt
    echo
done