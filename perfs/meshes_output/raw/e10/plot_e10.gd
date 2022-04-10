# Meta
set title "Dump every 10, * processes"

# Export
set terminal push
set terminal png size 1920,1080 enhanced font "Helvetica,12"
set output 'plot_e10.png'

# Settings
set grid
set autoscale 
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xlabel "Iterations"
set ylabel "Time (seconds)"

# Get files
list=system('ls -1 *e10_*.txt | tr "\n" " "')

# Plots
plot[2:] for [file in list] file using 1:2 with lines title file

set terminal pop
set output
replot
pause -1 ""