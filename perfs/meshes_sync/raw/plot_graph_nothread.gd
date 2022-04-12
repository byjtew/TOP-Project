# Meta
set title "Mesh synchronization: 1000 iterations, {20, 200} processes"

# Export
set terminal push
set terminal png size 1920,1080 enhanced font "Helvetica,12"
set output 'mesh_graph_nothread.png'

# Settings
set grid
set autoscale 
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set xlabel "Iterations"
set ylabel "Time (seconds)"

# Get files
graph_nothread_files=system('ls -1 timers-Graph-nothread*.txt | tr "\n" " "')

# Plots
plot [200:900] for [file in graph_nothread_files] file using 1:2 with lines title file

set terminal pop
set output
replot
pause -1 ""
