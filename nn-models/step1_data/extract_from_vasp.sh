#!/bin/sh

dir=$( head -1 "dir" | awk '(NR==1){print $1}' )
echo $dir
echo " "

DATA=`cat ./input_tag.dat`
while read tag
do
  echo $tag
  echo "---------------------------------"
  echo "Extracting data from VASP outputs"
  echo "---------------------------------"
  echo "DFT-MD data:" $tag
  echo " "


# BORN CHARGE FROM OUTCAR
  natom=$( grep NIONS ./data_vasp/$dir/OUTCAR_$tag | awk '{print $12}' )
#grepでNIONSが含まれる行を検索、awkでスペース区切りで左から12番目の文字列を抽出
#echo "natom" $natom

  echo "--" \
  >  ./data_nnp/born_$tag\.dat
  natom1=$(( $natom + 1 ))
  grep -A$natom1 "DIRECTION 3 BORN" ./data_vasp/$dir/OUTCAR_$tag \
  >> ./data_nnp/born_$tag\.dat
#DIRECTIONのある行からnatom1行後ろまでを出力。座標とborn_zの組み合わせ

# ENERGY FROM OSZICAR
  grep "F=" ./data_vasp/$dir/OSZICAR_$tag |\
  cat |\
  awk '{print substr($0, index($0, "F"), index($0, "E0") -1 )}' |\
  cut -d " " -f 2 \
  > ./data_nnp/energies_$tag\.dat


# POSITION FROM XDATCAR →そのまま
  cp ./data_vasp/$dir/XDATCAR_$tag ./data_nnp/positions_$tag\.dat


# Extract NUMBER OF IONS
  echo " " > ./data_nnp/info_$tag\.dat

  elem_a=$(  head -6 ./data_vasp/$dir/XDATCAR_$tag | tail -1 | awk '{print $1}' )
  natom_a=$( head -7 ./data_vasp/$dir/XDATCAR_$tag | tail -1 | awk '{print $1}' )

  elem_b=$(  head -6 ./data_vasp/$dir/XDATCAR_$tag | tail -1 | awk '{print $2}' )
  natom_b=$( head -7 ./data_vasp/$dir/XDATCAR_$tag | tail -1 | awk '{print $2}' )
  if [ -z "$natom_b" ]; then
    natom_b=0
  fi
  if [ -z "$elem_b" ]; then
    elem_b=none
  fi

  elem_c=$(  head -6 ./data_vasp/$dir/XDATCAR_$tag | tail -1 | awk '{print $3}' )
  natom_c=$( head -7 ./data_vasp/$dir/XDATCAR_$tag | tail -1 | awk '{print $3}' )
  if [ -z "$natom_c" ]; then
    natom_c=0
  fi
  if [ -z "$elem_c" ];then
    elem_c=none
  fi

  elem_d=$(  head -6 ./data_vasp/$dir/XDATCAR_$tag | tail -1 | awk '{print $4}' )
  natom_d=$( head -7 ./data_vasp/$dir/XDATCAR_$tag | tail -1 | awk '{print $4}' )
  if [ -z "$natom_d" ]; then
    natom_d=0
  fi
  if [ -z "$elem_d" ]; then
    elem_d=none
  fi

  natom_sum=$(( $natom_a + $natom_b + $natom_c + $natom_d ))
  echo "Total_number_of_ions" $natom_sum >> ./data_nnp/info_$tag\.dat
  echo " " >> ./data_nnp/info_$tag\.dat

  element_sum=0
  for n in a b c d  #`seq 1 $nelement`
  do
    dummy=$( eval echo '$'natom_$n )
    if [ $dummy -ne 0 ]; then
      element_sum=$(( $element_sum + 1 ))
    fi
  done

  echo "Total_ions:" $natom_sum
  echo "number_of_elements:" $element_sum

  #echo "Max number of elements ?"
  #read str
  #echo "Max_number_of_elements" $str \
  #>> ./data_nnp/info_$tag\_$temp\k.dat

  echo "Ions_per_type" $element_sum >> ./data_nnp/info_$tag\.dat
  echo " " >> ./data_nnp/info_$tag\.dat

  if [ $element_sum -ge 1 ]; then
    echo $elem_a $natom_a \
    >> ./data_nnp/info_$tag\.dat
    echo $elem_a $natom_a
  fi
  if [ $element_sum -ge 2 ]; then
    echo $elem_b $natom_b \
    >> ./data_nnp/info_$tag\.dat
    echo $elem_b $natom_b
  fi
  if [ $element_sum -ge 3 ]; then
    echo $elem_c $natom_c \
    >> ./data_nnp/info_$tag\.dat
    echo $elem_c $natom_c
  fi
  if [ $element_sum -ge 4 ]; then
    echo $elem_d $natom_d \
    >> ./data_nnp/info_$tag\.dat
    echo $elem_d $natom_d
  fi

  echo " " >> ./data_nnp/info_$tag\.dat

  steps=$( grep Direct ./data_vasp/$dir/XDATCAR_$tag | cat -n | tail -1 | awk '{print $1}' )
  echo "MD_steps" $steps

  echo "MD_steps" $steps >> ./data_nnp/info_$tag\.dat
  echo " " >> ./data_nnp/info_$tag\.dat


  echo "------------------------"
  echo "Finished extracting data"
  echo "------------------------"


done << FILE
$DATA
FILE
