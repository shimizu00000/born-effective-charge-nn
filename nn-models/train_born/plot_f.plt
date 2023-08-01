
unset key


plot "force_training.dat" u 1:4 pt 6 lt 1,\
     "force_training.dat" u 2:5 pt 6 lt 1,\
     "force_training.dat" u 3:6 pt 6 lt 1,\
     "force_test.dat" u 1:4 pt 8 lt 3,\
     "force_test.dat" u 2:5 pt 8 lt 3,\
     "force_test.dat" u 3:6 pt 8 lt 3,\
     x
