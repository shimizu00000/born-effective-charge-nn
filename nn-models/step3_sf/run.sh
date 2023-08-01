#!/bin/sh
#PBS -l nodes=4:ppn=24
#PBS -N sf

cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
date > testjob.dat
echo 'PBS_O_WORKDIR: '$PBS_O_WORKDIR >> testjob.dat
echo 'PBS_NODEFILE: '$PBS_NODEFILE >> testjob.dat
echo 'PBS_JOBID: '$PBS_JOBID >> testjob.dat
pwd >> testjob.dat

mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 24 ./a.out >> testjob.dat

date >> testjob.dat
