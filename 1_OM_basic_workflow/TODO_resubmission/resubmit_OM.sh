#!/bin/bash
#SBATCH --job-name=rerun_unfinished_jobs
#SBATCH --account=smith
#SBATCH -o /scicore/home/smith/GROUP/smc_lai/JOB_OUT/resubmit_OM.out
#SBATCH --mem=2G
#SBATCH --qos=30min
#SBATCH --cpus-per-task=1
###########################################
# Script for resubmitting OpenMalaria simulations as jobs to the cluster
# Arguments:
#		OM_DIR: directory containing all the OpenMalaria simulation output files
#		SCENARIO_FOLDER: directory with all the OpenMalaria scenario .xml files
# Calling the script:
#	sbatch --array=1-NUM resubmit_OM.sh OM_DIR SCENARIO_FOLDER
#
# created on 01.02.2019
# monica.golumbeanu@unibas.ch
###########################################

OM_DIR=$1
SCENARIO_FOLDER=$2

# Load OpenMalaria module and change to folder with resource files
module purge
ml OpenMalaria/38.0-goolf-1.7.20-Python-2.7.11
cd /scicore/home/smith/GROUP/smc_lai/OM_schema38

# module load OpenMalaria/32-RC3-goolf-1.7.20-Python-2.7.11
# cd /scicore/home/smith/golmon00/OM_schema32

# Select scenario file in array
scenario_files=(${SCENARIO_FOLDER}*.xml)
ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)
scenario_file=${scenario_files[$ID]}
#echo "Running simulation for $scenario_file"

# Check if the simulation outputs exist and if not, rerun OpenMalaria
scenario_name=$(basename "$scenario_file" ".xml")
#echo $scenario_name
output1=$scenario_name"_out.txt"
file_count=$(find $OM_DIR -name $output1 | wc -l)

if [[ $file_count -eq 0 ]]; then
    echo "Running simulation for $scenario_file"
    output1=$OM_DIR$scenario_name"_out.txt"
    output2=$OM_DIR$scenario_name"_cts.txt"
    echo "Outputs will be written to $output1 and $output2"
    openMalaria --scenario $scenario_file --output $output1 --ctsout $output2
    echo "OpenMalaria simulation ended."
    parentdir="$(dirname "$OM_DIR")"
    echo "simulation ended for $scenario_name $ID" >> $parentdir"/resubmit.txt"
fi
