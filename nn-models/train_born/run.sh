#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -N jobname

cd $PBS_O_WORKDIR
ulimit -s unlimited
NPROCS=`wc -l <$PBS_NODEFILE`
date >> testjob.dat
hostname >> testjob.dat
pwd >> testjob.dat
echo 'PBS_O_WORKDIR: '$PBS_O_WORKDIR >> testjob.dat
echo 'PBS_NODEFILE: '$PBS_NODEFILE >> testjob.dat
echo 'PBS_JOBID: '$PBS_JOBID >> testjob.dat

mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 16 ./a.out >> output.dat

date >> testjob.dat
