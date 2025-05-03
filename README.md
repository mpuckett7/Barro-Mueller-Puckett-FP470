# CS 470 Final Project

Authors: Gavin Barro, Mason Puckett, Beau Mueller

Running Experiments:

NOTE:  run make clean in the examples's (rayleighBenard and bifurcation3d) directory and OpenLB folder and run module load mpi

MPI (SBATCH) ->
- For MPI use look under /MPI-TestingScripts and under whichever example you wish to run. Then give yourself execute permission for the launch script and run the launch script only. If you wish to cat all of the results to your terminal feel free to use the view script however that is A LOT of text and I would recommend just using the text files themselves. 
- Sometimes the sbatch jobs can take an absurdly long time to complete because when running with only MPI a process can slip away so if you see anything running for over like 5 minutes go ahead and cancel it. It will sometimes still give you all the information about the simulation except the mpimanager missed wrapping a process up after the simulation ended

OMP (NOT SBATCH) ->
- For OMP use look under /OMP-TestingScripts and under whichever example you wish to run. Then give yourself execute permission for the run script and execute the run script. This one will print out the make info onto your terinmal and then you will have to wait for all the tests to be done before getting your terminal back. (./SCRIPTNAME to run)
- Output for all the tests can be found in text files in the examples directory under OpenLB. For example, /OpenLB/examples/thermal/raylegihBenard2d is the path to where the output files are from the top of the repository. 

Hybrid (SBATCH) -> 
- The hybrid testing scripts are set up in the same manner the MPI testing scripts are. Give yourself execute permission to the launch script and then ./ the launch script

Weak Scaling (NOT SBATCH) ->
- To run our version of weak scaling, tuning the model resolution, then you can find a script under /Weak-Scaling. We only have a weak scaling script for running the experiments, simply pick the example and run the script. (./SCRIPTNAME)
- Output for thise script has been sent to text files within the example directory as well, you will not be able to use your terminal while this is executing either.

Contribution:
- We have decided to add to this repository by writing a .hh file that allows for portable code adjustment of the model resolution. It is under /src/utilities as resolutionInput.hh 
- As stated in the weak scaling experiments some examples have the option to use model resolution as a command line parameter already while others don't so we went ahead and created a method that you can call from any example which allows you to pass in the model resolution as a command line argument using a -r flag.

