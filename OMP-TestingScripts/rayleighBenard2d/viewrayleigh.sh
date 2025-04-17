#!/bin/bash
# Run to view full or partial results for OpenMP testing of the rayleighBenard2d example

for n in 1 2 4 8 16 32 64; do
    echo "== rayleighBenard2d with $n threads =="
    cat rayleigh_changed-$n.txt
    echo
done