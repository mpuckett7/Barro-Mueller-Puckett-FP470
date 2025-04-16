#!/bin/bash
# Run to view full or partial results for MPI testing of the rayleighBenard2d example


# for n in 1 2  4 8 16 32 64; do
#     echo "== rayleighBenard2d with $n processes, NO changes =="
#     cat rayleigh_unchanged-$n.txt
#     echo
# done

for n in 1 2 4 8 16 32 64; do
    echo "== rayleighBenard2d with $n processes, with contributions =="
    cat rayleigh_changed-$n.txt
    echo
done