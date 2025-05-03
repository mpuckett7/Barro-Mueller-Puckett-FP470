#!/bin/bash

for n in 2 4 8; do
    echo "== Bifurcation, LaGrange with $n processes and 32 threads/process =="
    cat bifurcation3d-hybrid-$n.txt
    echo
done