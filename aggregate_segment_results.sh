#!/bin/sh

# Directives
#PBS -N aggregate_segment_results
#PBS -W group_list=yetistats
#PBS -l nodes=1,walltime=2:00:00,mem=6000mb
#PBS -M mo2499@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories (SSCC example here)
#PBS -o localhost:/vega/stats/users/mo2499/bbSimulation/outputLog/
#PBS -e localhost:/vega/stats/users/mo2499/bbSimulation/errorLog/

# Run MATLAB function
matlab -nosplash -nodisplay -nodesktop -r "aggregate_segment_results" > output5

# End of script
