#!/bin/sh

# Directives
#PBS -N segmentImage
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=22:00:00,mem=10000mb
#PBS -t 1-12,17-24,29-36,41-48
#PBS -M mo2499@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories (SSCC example here)
#PBS -o localhost:/vega/stats/users/mo2499/bbSimulation/outputLog/
#PBS -e localhost:/vega/stats/users/mo2499/bbSimulation/errorLog/

# Run MATLAB function
matlab -nosplash -nodisplay -nodesktop -r "simulate_par_segmentImage($PBS_ARRAYID)" > output4

# End of script
