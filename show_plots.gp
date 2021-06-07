#######################
# Plotting CFA curves #
#######################
set terminal wxt 0
set title "CFA - Fitted"
plot "cfa_spectra.dat" u 1:2 lw 2 lc 7 w l title "CFA Red", \
"cfa_spectra.dat" u 1:3 lw 2 lc 10 w l title "CFA Green", \
"cfa_spectra.dat" u 1:4 lw 2 lc 14 w l title "CFA Blue", \
"cfa_ref.dat" u 1:2 axis x1y2 lc 7 dt 4 w l title "Reference Red", \
"cfa_ref.dat" u 1:3 axis x1y2 lc 10 dt 4 w l title "Reference Green", \
"cfa_ref.dat" u 1:4 axis x1y2 lc 14 dt 4 w l title "Reference Blue"
reset

##################
# Plotting costs #
##################
set terminal wxt 1
set title "Cost"
plot "cfa_fitting.dat" u 1:3 lc 7 dt 1 w l title "Initial Cost", \
"cfa_fitting.dat" u 1:4 lc 10 dt 1 w l title "Final Cost"
reset

###################
# Plotting Errors #
###################
set terminal wxt 2
set multiplot

set lmargin at screen 0.075
set rmargin at screen 0.975
set bmargin at screen 0.1
set tmargin at screen 0.95
set size ratio 1.125
set xrange [0:979]
set yrange [0:1102]
unset tics
plot "CIExy1931_bg.png" binary filetype=png w rgbimage notitle

set size ratio 1.125
set tics
set xrange [0:0.8]
set yrange [0:0.9]
plot "cfa_vector.dat" u 1:2:($3-$1):($4-$2) w vec head filled lt 2 lc 0 title "Matrix Errors"

unset multiplot
reset

####################
# Plotting Inverse #
####################
set terminal wxt 3
set multiplot

set lmargin at screen 0.075
set rmargin at screen 0.975
set bmargin at screen 0.1
set tmargin at screen 0.95
set size ratio 1.125
set xrange [0:979]
set yrange [0:1102]
unset tics
plot "CIExy1931_bg.png" binary filetype=png w rgbimage notitle

set size ratio 1.125
set tics
set xrange [0:0.8]
set yrange [0:0.9]
plot "inverse.dat" u 3:4:($1-$3):($2-$4) w vec head filled lt 2 lc 0 title "Inverse"

unset multiplot
reset
