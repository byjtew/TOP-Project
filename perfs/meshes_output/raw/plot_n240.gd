# Meta
set title "Dump every 100, * processes"

# Export
set terminal push
set terminal png size 1920,1080 enhanced font "Helvetica,12"
set output 'plot_n240.png'

# Settings
set grid
set autoscale 
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xlabel "Iterations"
set ylabel "Time (seconds)"

# Get files
recv_list=system('ls -1 Recv*n240*.txt | tr "\n" " "')
irecv_list=system('ls -1 Irecv*n240*.txt | tr "\n" " "')
gather_list=system('ls -1 Gather*n240*.txt | tr "\n" " "')

# Plots
plot  [2:] for [file in recv_list] file using 1:2 with lines lt rgb "red" title file, \
for [file in irecv_list] file using 1:2 with lines lt rgb "blue" title file, \
for [file in gather_list] file using 1:2 with lines lt rgb "goldenrod" title file


set terminal pop
set output
replot
pause -1 ""
