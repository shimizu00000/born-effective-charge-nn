#!/bin/sh 
#PBS -l nodes=8:ppn=24
#PBS -N testjob
cd $PBS_O_WORKDIR
ulimit -s unlimited
NPROCS=`wc -l <$PBS_NODEFILE`
date > testjob.dat
hostname >> testjob.dat
pwd >> testjob.dat
echo 'PBS_O_WORKDIR: '$PBS_O_WORKDIR >> testjob.dat
echo 'PBS_NODEFILE: '$PBS_NODEFILE >> testjob.dat
echo 'PBS_JOBID: '$PBS_JOBID >> testjob.dat

mpiexec -iface ib0 -launcher rsh -machinefile $PBS_NODEFILE -ppn 24 /home/shimizu/lammps-12Dec18/build-hara/lmp -partition 8x24 -in in.neb >> testjob.dat
