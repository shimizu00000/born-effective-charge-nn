#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -N jobname

cd $PBS_O_WORKDIR
ulimit -s unlimited
NPROCS=`wc -l <$PBS_NODEFILE`
date >> output.dat
hostname >> output.dat
pwd >> output.dat
echo 'PBS_O_WORKDIR: '$PBS_O_WORKDIR >> output.dat
echo 'PBS_NODEFILE: '$PBS_NODEFILE >> output.dat
echo 'PBS_JOBID: '$PBS_JOBID >> output.dat

mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 ./a.out >> output.dat
