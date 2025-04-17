#!/bin/bash
# Run to view full or partial results for MPI testing of the rayleighBenard2d example

    echo "== rayleighBenard2d with 1 processes =="
    cat rayleigh1-1.txt
    echo

    echo "== rayleighBenard2d with 2 processes =="
    cat rayleigh2-2.txt
    echo

for n in 4 8 16 32; do
    echo "== rayleighBenard2d with $n processes =="
    cat rayleigh_changed-$n.txt
    echo
done