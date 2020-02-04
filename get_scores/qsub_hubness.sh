#!/bin/bash

## Torque Configuration
# ==============================================================================

# resources
#PBS -l walltime=150:00:00
#PBS -l mem=20gb
#PBS -l nodes=1:ppn=1
#PBS -q batch
 
# Information
#PBS -N hubness_satija_data
#PBS -M elise.amblard@curie.fr
#PBS -m abe
 
# Other
#PBS -j oe
#PBS -o /data/tmp/zela/ML/hubness
#PBS -p 900


## Command section
# ==============================================================================

cd $PBS_O_WORKDIR
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

# Main
# ==============================================================================
export PATH="/data/users/eamblard/tmp/R-3.5.3/bin:$PATH"
script_dir="/data/tmp/zela/ML/hubness"
script="1_get_hubness_scores_minkow.R"
working_directory="/data/tmp/zela/ML/hubness"

echo Start = `date`
cd $working_directory
Rscript ${script_dir}/${script}

echo End = `date`