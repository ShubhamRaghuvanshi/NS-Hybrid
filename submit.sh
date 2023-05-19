#!/bin/bash
# Submission script: "tasks are evenly distributed across nodes"

# Job name
#SBATCH --job-name=NS_hybrid_test
# Output file name
#SBATCH --output=test.out
#SBATCH --error=test.err
#
# Set the required partition [change]
#SBATCH --partition=parallel-short

# Number of processes
#SBATCH --ntasks=8

# Process distribution per node
#SBATCH --ntasks-per-node=2

# Number of nodes
#SBATCH --nodes=4

##SBATCH --sockets-per-node=1

#SBATCH --cpus-per-task=16

# Memory per process
#SBATCH --mem-per-cpu=4000
#
# Total wall-time
#SBATCH --time=08:05:00
#
# The below statement is required if the code is floating-point intensive and CPU-bound [Optional]
#SBATCH --threads-per-core=1
#


SimFile=SimulationDetails.log
if test -f "$SimFile"; then  
rm $SimFile  
fi

echo 'SLURM_NPROCS : ' $SLURM_NPROCS >> $SimFile
echo 'CPU per process :' $SLURM_CPUS_PER_TASK >> $SimFile

curdir=$(pwd)
echo 'IO Directory :' $curdir >> $SimFile

start_time=$SECONDS

mpirun --mca btl_openib_allow_ib 1 --bind-to none -x OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK $curdir/part.exe $curdir > $curdir/simulation.log
#mpirun -np $SLURM_NPROCS --bind-to none -x OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK $curdir/part.exe $curdir > $curdir/simulation.log

end_time=$SECONDS
elapsed=$(( end_time - start_time )) 

echo 'Time of execution : '$elapsed >> SimulationDetails.log







