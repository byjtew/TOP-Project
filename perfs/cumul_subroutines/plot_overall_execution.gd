# Meta
set title "Cumulated time in subroutines: 10.000 iterations, 100 processes, 1 IO every 100 iterations"

# Style
set style data histogram
set style fill solid border
set style histogram clustered

# Export
set terminal push
set terminal png size 1920,1080 enhanced font "Helvetica,12"
set output 'cumul_subroutines_histogram.png'


# Settings
set ylabel "Time (seconds)"

# Plots
plot  'timers-cumul-subroutine.txt' using 2:xticlabels(1) title columnheader

set terminal pop
set output
replot
pause -1 ""
