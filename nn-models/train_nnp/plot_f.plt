set terminal png
set output "f_train.png"

set size 0.9
set size square

#set xrange [-12:12]
#set yrange [-12:12]

set xtics 10 
set ytics 10


unset key


plot "force_training.dat" u 1:4 pt 6 lt 1,\
     "force_training.dat" u 2:5 pt 6 lt 1,\
     "force_training.dat" u 3:6 pt 6 lt 1,\
     x
