#######################
# Plotting CFA curves #
#######################
set terminal pngcairo
set output "cfa_fitted.png"
plot "cfa_spectra.dat" u 1:2 lw 2 lc 7 w l title "Rot", \
"cfa_spectra.dat" u 1:3 lw 2 lc 10 w l title "Grün", \
"cfa_spectra.dat" u 1:4 lw 2 lc 14 w l title "Blau"
reset

############################
# Plotting Iterations used #
############################
set terminal pngcairo
set output "cfa_iterations.png"
plot "cfa_fitting.dat" u 1:2 lw 2 lc 0 dt 1 w l title "Iterationen"
reset

##################
# Plotting costs #
##################
set terminal pngcairo
set output "cfa_cost.png"
plot "cfa_fitting.dat" u 1:3 lc 7 dt 1 w l title "Initiale Kosten", \
"cfa_fitting.dat" u 1:4 lc 10 dt 1 w l title "Finale Kosten"
reset

###################
# Plotting Errors #
###################
set terminal pngcairo size 640,720
set output "cfa_errors.png"
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
plot "cfa_vector.dat" u 1:2:($3-$1):($4-$2) w vec head filled lt 2 lc 0 title "Matrixfehler"

unset multiplot
reset

#####################
# Plotting Sampling #
#####################
set terminal pngcairo size 640,720
set output "cfa_sampling.png"
set multiplot

set lmargin at screen 0.075
set rmargin at screen 0.975
set bmargin at screen 0.1
set tmargin at screen 0.95
set xrange [0:979]
set yrange [0:1102]
unset tics
plot "CIExy1931_bg.png" binary filetype=png w rgbimage

set tics
set xrange [0:0.8]
set yrange [0:0.9]
plot "cfa_sampling.dat" title "Samples"

unset multiplot
reset

####################
# Plotting Inverse #
####################
set terminal pngcairo size 640,720
set output "cfa_inverse.png"
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
