#!/bin/bash
# Run to view full or partial results for OpenMP testing of the EXAMPLE-NAME example


for n in 1 2  4 8 16 32 64; do
    echo "== EXAMPLE-NAME with $n threads, NO changes =="
    cat EXAMPLE-NAME_unchanged-$n.txt
    echo
done

for n in 1 2 4 8 16 32 64; do
    echo "== EXAMPLE-NAME with $n threads, with contributions =="
    cat EXAMPLE-NAME_changed-$n.txt
    echo
done