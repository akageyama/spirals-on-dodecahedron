#
#



set xrange [0:2*pi]
set yrange [0:pi]

set isosamples 400,200

set view 0, 0, 1

splot 'output.data' with pm3d

pause -1

# set terminal png
set terminal png size 1280, 960
set output 'output.png'
splot 'output.data' with pm3d
#replot
