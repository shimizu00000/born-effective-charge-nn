set terminal png
set output "e_train.png"

set size 0.9 
set size square

#set xrange [-6.3:-5.5]
#set yrange [-6.3:-5.5]
set xtics 0.2 
set ytics 0.2

unset key


plot "energy_training.dat" u 1:2 pt 6 lt 1,\
     x 
