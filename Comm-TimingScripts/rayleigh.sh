cd ../OpenLB
make clean
module load mpi
make

cd examples/thermal/rayleighBenard2d
make clean && make

echo "testing 1-1"
salloc -Qn 1 --nodes=1 --ntasks-per-node=1 mpirun ./rayleighBenard2d > comtime-1-1.out
echo "testing 1-2"
salloc -Qn 2 --nodes=1 --ntasks-per-node=2 mpirun ./rayleighBenard2d > comtime-1-2.out
echo "testing 1-4"
salloc -Qn 4 --nodes=1 --ntasks-per-node=4 mpirun ./rayleighBenard2d > comtime-1-4.out
echo "testing 2-1"
salloc -Qn 2 --nodes=2 --ntasks-per-node=1 mpirun ./rayleighBenard2d > comtime-2-1.out
echo "testing 2-2"
salloc -Qn 4 --nodes=2 --ntasks-per-node=2 mpirun ./rayleighBenard2d > comtime-2-2.out
echo "testing 2-4"
salloc -Qn 8 --nodes=2 --ntasks-per-node=4 mpirun ./rayleighBenard2d > comtime-2-4.out
echo "testing 4-1"
salloc -Qn 4 --nodes=4 --ntasks-per-node=1 mpirun ./rayleighBenard2d > comtime-4-1.out
echo "testing 4-2"
salloc -Qn 8 --nodes=4 --ntasks-per-node=2 mpirun ./rayleighBenard2d > comtime-4-2.out
echo "testing 4-4"
salloc -Qn 16 --nodes=4 --ntasks-per-node=4 mpirun ./rayleighBenard2d > comtime-4-4.out
echo "testing 8-1"
salloc -Qn 8 --nodes=8 --ntasks-per-node=1 mpirun ./rayleighBenard2d > comtime-8-1.out
echo "testing 8-2"
salloc -Qn 16 --nodes=8 --ntasks-per-node=2 mpirun ./rayleighBenard2d > comtime-8-2.out
echo "testing 8-4"
salloc -Qn 32 --nodes=8 --ntasks-per-node=4 mpirun ./rayleighBenard2d > comtime-8-4.out
