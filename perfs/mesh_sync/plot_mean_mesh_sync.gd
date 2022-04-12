# Meta
set title "Mesh synchronization: 30.000 iterations, 100 processes"

# Style
set style data histogram
set style fill solid border
set style histogram clustered

# Export
set terminal push
set terminal png size 1920,1080 enhanced font "Helvetica,12"
set output 'mesh_sync_histogram.png'


# Settings
set ylabel "Time (seconds)"

# Plots
plot for [COL=2:4] 'timers-mean-MeshSync.txt' using COL:xticlabels(1) title columnheader

set terminal pop
set output
replot
pause -1 ""
