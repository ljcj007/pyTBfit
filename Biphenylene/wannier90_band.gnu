set data style dots
set nokey
set xrange [0: 4.16328]
set yrange [-23.60707 :  7.40298]
set arrow from  0.83843, -23.60707 to  0.83843,   7.40298 nohead
set arrow from  1.53623, -23.60707 to  1.53623,   7.40298 nohead
set arrow from  2.62705, -23.60707 to  2.62705,   7.40298 nohead
set arrow from  3.32486, -23.60707 to  3.32486,   7.40298 nohead
set xtics (" G "  0.00000," X "  0.83843," M "  1.53623,"G/X"  2.62705," Y "  3.32486," M "  4.16328)
 plot "wannier90_band.dat"
