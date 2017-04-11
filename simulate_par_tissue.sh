#!/bin/sh

# Directives
#PBS -N simulate_tissue
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=20:00:00,mem=20000mb
#PBS -t 1-12
#PBS -M mo2499@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories (SSCC example here)
#PBS -o localhost:/vega/stats/users/mo2499/bbSimulation/outputLog/
#PBS -e localhost:/vega/stats/users/mo2499/bbSimulation/errorLog/

# Run MATLAB function
matlab -nosplash -nodisplay -nodesktop -r "simulate_par_tissue($PBS_ARRAYID)" > output0

# End of script
