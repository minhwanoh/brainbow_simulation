#!/bin/sh

# Directives
#PBS -N segmentImageWithIncorrectCellCount
#PBS -W group_list=yetistats
#PBS -l nodes=1:ppn=1,walltime=45:00:00,mem=6000mb
#PBS -M us2157@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories (SSCC example here)
#PBS -o localhost:/vega/stats/users/us2157/bb/outputLog/
#PBS -e localhost:/vega/stats/users/us2157/bb/errorLog/

# Run MATLAB function
matlab -nosplash -nodisplay -nodesktop -r "simulate_segmentImageWithIncorrectCellCount" > output1

# End of script
