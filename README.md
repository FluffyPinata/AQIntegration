# AQIntegration

To run each file, first compile them with the following commands:
$ g++ -O -fopenmp aq-sequential.cpp -o aq-seq
$ g++ -O -fopenmp aq-openmp.cpp -o aq-omp
$ mpicxx aq-mpi.cpp -o aq-mpi 

Then initialize a job with:
$  sbatch c_slurm.sh
$  sbatch mpi_slurm.sh
$  sbatch openmp_slurm.sh

It may fail on the initial run after compliation, if it does just run it again and it will work.


