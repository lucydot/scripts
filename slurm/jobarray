#!/bin/bash --login
#PBS -N 142
#PBS -l select=32
#PBS -l walltime=00:03:00
#PBS -A e05-gener-wal

FOLDERS=(0.08 0.04 0.02 0.01 0.0075 0.005 0.0025 0.001)
JOBS=/work/e05/e05/lucy/jobs/142
MPI=96

# Note: select = No. Nodes (24 CPUS per node) # Also from the architecture the sweet spots should be:
# 4 nodes (96 cores), which is one blade
# and multiples of four up to
# 64 nodes (1536 cores), which is one chassis.

# Hey lucy, if you need to c .. make sure you c - instead

NFOLDERS=${#FOLDERS[@]}

# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

# Change to the directory that the job was submitted from.
cd $PBS_O_WORKDIR

# Set the number of threads to 1
# This prevents any system libraries from automatically 
# using threading.
export OMP_NUM_THREADS=1

# VASP makes big files
ulimit -s unlimited

# Number of MPI processes
NPROCT=`qstat -f $PBS_JOBID | awk '/Resource_List.mpiprocs/ {print $3}'`

for f in ${FOLDERS[@]}
do
    cd $JOBS/$f
    echo "Submitting job $JOBS/$FOLDERS dotty...."    
    aprun -n $MPI /home/e05/e05/lucy/bin/vasp_std > vasp.out &
done

wait


