#!/bin/bash
echo "LAMMPS dump-file name? (preset)"
filename="dump.melt"
echo "Reading $filename ..."


natom=$( sed -n 4p $filename | awk '(NR==1){print $1}' )
echo "Number of atoms = $natom"

nstep=$( grep TIMESTEP $filename | wc | awk '(NR==1){print $1}' )
echo "Number of MD steps = $nstep"

line=$(( $natom * 3 ))
head -"$line" $filename | grep -A1 TIMESTEP > ./dummy
step0=$( sed -n 2p ./dummy | awk '(NR==1){print $1}' )
step1=$( sed -n 5p ./dummy | awk '(NR==1){print $1}' )
dgap=$(( $step1 - $step0 ))
tail -"$line" $filename | grep -A1 TIMESTEP > ./dummy
stepf=$( tail -1 ./dummy | awk '(NR==1){print $1}' )
echo "initial: $step0, next: $step1, last: $stepf (every $dgap steps)"
rm ./dummy


dline=$(( $natom + 9 ))
for n in `seq $nstep`; do

  line1=$(( 1 + ($dline*($n-1)) ))
  line2=$(( $line1 + $natom + 8 ))
  sed -n "$line1","$line2"p $filename > ./dummy
  echo $n"/"$nstep  $line1 $line2    #1~585, 586~1170, #576+9,586

  cstep=$( sed -n 2p ./dummy | awk '(NR==1){print $1}' )

  echo "HEADER $cstep" > ./data.relax
  echo "$natom atoms" >> ./data.relax
  echo "3 atom types" >> ./data.relax

  xlo=$( sed -n 6p ./dummy | awk '(NR==1){print $1}' ); xhi=$( sed -n 6p ./dummy | awk '(NR==1){print $2}' )
  ylo=$( sed -n 7p ./dummy | awk '(NR==1){print $1}' ); yhi=$( sed -n 7p ./dummy | awk '(NR==1){print $2}' )
  zlo=$( sed -n 8p ./dummy | awk '(NR==1){print $1}' ); zhi=$( sed -n 8p ./dummy | awk '(NR==1){print $2}' )

  echo $xlo $xhi "xlo xhi" >> ./data.relax
  echo $ylo $yhi "ylo yhi" >> ./data.relax
  echo $zlo $zhi "zlo zhi" >> ./data.relax
  echo ""           >> ./data.relax
  echo "Masses"     >> ./data.relax
  echo ""           >> ./data.relax
  echo "1 7.010 Li" >> ./data.relax
  echo "2 30.974 P" >> ./data.relax
  echo "3 32.066 S" >> ./data.relax
  echo ""           >> ./data.relax
  echo "Atoms"      >> ./data.relax
  echo ""           >> ./data.relax 

  xlof=`printf "%f" $xlo`; xhif=`printf "%f" $xhi`
  ylof=`printf "%f" $ylo`; yhif=`printf "%f" $yhi`
  zlof=`printf "%f" $zlo`; zhif=`printf "%f" $zhi`

  xlat=`echo "scale=20; $xhif - $xlof" | bc -l`
  ylat=`echo "scale=20; $yhif - $ylof" | bc -l`
  zlat=`echo "scale=20; $zhif - $zlof" | bc -l`

  for nn in `seq $natom`; do
    latom=$(( 9 + $nn ))
    #echo $latom
    sed -n "$latom"p ./dummy > ./dummya
    atomID=$( tail -1 ./dummya | awk '(NR==1){print $1}' )
    elemID=$( tail -1 ./dummya | awk '(NR==1){print $2}' )
    posisx=$( tail -1 ./dummya | awk '(NR==1){print $3}' ); posisxf=`printf "%f" $posisx`
    posisy=$( tail -1 ./dummya | awk '(NR==1){print $4}' ); posisyf=`printf "%f" $posisy`
    posisz=$( tail -1 ./dummya | awk '(NR==1){print $5}' ); posiszf=`printf "%f" $posisz`
    
    posix=`echo "scale=20; $xlat * $posisxf" | bc -l`
    posiy=`echo "scale=20; $ylat * $posisyf" | bc -l`
    posiz=`echo "scale=20; $zlat * $posiszf" | bc -l`

    echo $atomID $elemID $posix $posiy $posiz >> ./data.relax

  done
  rm ./dummya

  echo "Execute lmp2pos.py ..."
  python ./lmp2pos.py
  sed -e "s/H   He  Li/Li P S/g" ./POSCAR.relax > ./POSCAR."$cstep"
  echo "Finish!"

done
rm ./dummy
