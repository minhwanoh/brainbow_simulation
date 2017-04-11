#!/bin/sh

# Directives
#PBS -N segmentImage
#PBS -W group_list=yetistats
#PBS -l nodes=1:ppn=1,walltime=45:00:00,mem=6000mb
#PBS -M mo2499@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories (SSCC example here)
#PBS -o localhost:/vega/stats/users/mo2499/bbSimulation/outputLog/
#PBS -e localhost:/vega/stats/users/mo2499/bbSimulation/errorLog/

# Run MATLAB function
matlab -nosplash -nodisplay -nodesktop -r "simulate_segmentImage" > output1

# End of script
