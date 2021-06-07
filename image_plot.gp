#######################
# Plotting Chromacity #
#######################
set terminal pngcairo size 640,720
set output "image_spectra.png"
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
plot "image_spectra.dat" u 1:2 lc 16 pt 17 title "xy color"
