#/bin/bash

for n in 2 4 8 16 32; do
    echo "== RaleighBenard2d with $n processes and 32 threads/process =="
    cat rayleigh-hybrid-$n.txt
    echo
done
