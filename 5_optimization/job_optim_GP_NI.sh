#!/bin/bash
#SBATCH --job-name=optim_GP
#SBATCH --account=smith
#SBATCH -o /scicore/home/penny/burlyd00/smc_lai/JOB_OUT/optimisation_analysis.out
#SBATCH --mem=10G
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=1

#######################################
# script for optimisation analysis
#
# created 08.01.2020
# monica.golumbeanu@unibas.ch
######################################
ml purge
ml R/3.6.0-foss-2018b

GP_DIR1=$1
GP_DIR2=$2
PARAM_RANGES_FILE=$3
OPT_DEST_DIR=$4
OPT_SETUP_FILE=$5
ROW=$6

# IMPORTANT: the number of files must equal to the task array length (index starts at 0)
gp_files1=(${GP_DIR1}*.RData)
gp_files2=(${GP_DIR2}*.RData)

# Select scenario file in array
ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)
gp_file1=${gp_files1[$ID]}
gp_file2=${gp_files2[$ID]}

echo "Optimisation for $gp_file1"

Rscript get_optim_profile_noninferiority.R $gp_file1 $gp_file2 $PARAM_RANGES_FILE $OPT_DEST_DIR $OPT_SETUP_FILE $ROW

