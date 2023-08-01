#!/bin/bash
rm energy_training.dat energy_test.dat
rm force_training.dat  force_test.dat

for num_mpi in `seq 1 32`
do

 a=$(( num_mpi - 1  ))
 if [ $a -lt 10 ]; then
   cat energy_training000$a\.dat >> energy_training.dat
   cat energy_test000$a\.dat >> energy_test.dat
   cat force_training000$a\.dat >> force_training.dat
   cat force_test000$a\.dat >> force_test.dat
 elif [ $a -ge 10 -a $a -lt 100 ];then
   cat energy_training00$a\.dat >> energy_training.dat
   cat energy_test00$a\.dat >> energy_test.dat
   cat force_training00$a\.dat >> force_training.dat
   cat force_test00$a\.dat >> force_test.dat
 elif [ $a -ge 100 -a $a -lt 1000 ];then
   cat energy_training0$a\.dat >> energy_training.dat
   cat energy_test0$a\.dat >> energy_test.dat
   cat force_training0$a\.dat >> force_training.dat
   cat force_test0$a\.dat >> force_test.dat
 fi

done
