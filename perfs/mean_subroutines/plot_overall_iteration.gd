# Meta
set title "Mean time in subroutines: 30.000 iterations, 100 processes"

# Style
set style data histogram
set style fill solid border
set style histogram clustered

# Export
set terminal push
set terminal png size 1920,1080 enhanced font "Helvetica,12"
set output 'mean_subroutines_histogram.png'


# Settings
set ylabel "Time (seconds)"

# Plots
plot  'timers-mean-subroutine.txt' using 2:xticlabels(1) title columnheader

set terminal pop
set output
replot
pause -1 ""
