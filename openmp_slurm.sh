#!/bin/bash
##############################################################################
# FILE: openmp_slurm
# DESCRIPTION:
#   Simple Slurm job command script for openmp code on Heracles
# AUTHOR: Manh and Iris
# LAST REVISED: 08/2020
##############################################################################
### Job commands start here
## echo '=====================JOB STARTING=========================='
#SBATCH --partition=day-long-cpu        ### Partition
## PartitionName=day-long-cpu Priority=20000 Default=NO MaxTime=1-0:00:00 State=UP Nodes=node[2-17]
#SBATCH --job-name=openmpcode              ### Job Name
#SBATCH --output=slurm_output.%j        ### File in which to store job output
#SBATCH --error=slurm_error.%j          ### File in which to store job error messages
#SBATCH --time=0-00:01:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job, default = 1
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node, default = 1
#SBATCH --cpus-per-task=48  ### Number of threads per task (OMP threads) and set SLURM_CPUS_PER_TASK
#SBATCH --exclusive             ### no shared resources within a node
#SBATCH --mail-type=ALL                   ### email alert at start, end and abortion of execution
## SBATCH --mail-user=myemail ### send mail to this address (uncommnet this line in case using email)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo -n 'cpus per task : '; echo $SLURM_CPUS_PER_TASK
echo -n 'omp num threads : ';  echo $OMP_NUM_THREADS
## slurm will decide which node the code will run according to the partition ste in the job. No need ssh node on heracles
./mm-omp "$@"
## the script receives the matrix size from keyboard
## execute this script like that
## sbatch openmp_slurm.sh size  ###   size is an integer number for the matrix size
### Display some diagnostic information
echo ' '
echo '=====================JOB DIAGNOTICS========================'
date
echo -n 'Job is executed on ';hostname
echo -n 'My jobid is '; echo $SLURM_JOBID
#echo 'My path is:' 
#echo $PATH
sinfo -s
echo 'My job info:'
squeue -j $SLURM_JOBID
echo 'Machine info'
echo ' '
echo '========================ALL DONE==========================='
