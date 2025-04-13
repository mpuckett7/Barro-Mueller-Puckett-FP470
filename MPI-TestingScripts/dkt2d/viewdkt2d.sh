#!/bin/bash
# Run to view full or partial results for MPI testing of the dkt2d example


for n in 1 2  4 8 16 32 64; do
    echo "== DKT2D with $n processes, NO changes =="
    cat dkt2d_unchanged-$n.txt
    echo
done

for n in 1 2 4 8 16 32 64; do
    echo "== DKT2D with $n processes, with contributions =="
    cat dkt2d_changed-$n.txt
    echo
done