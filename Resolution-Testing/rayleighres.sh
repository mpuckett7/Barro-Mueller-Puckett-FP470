cd ../
cd OpenLB/
module load mpi
make clean
make
cd examples/thermal/rayleighBenard2d


for t in 21 22 23 24 25 26 27 28 29 30; do
    make clean
    make
    OMP_NUM_THREADS=12 salloc -Q -n 2 mpirun ./rayleighBenard2d -r $t | grep -A 5 '\[Timer\] ----------------Summary:Timer----------------'&> ../../../../Resolution-Testing/ray_res_r=$t.txt
done
