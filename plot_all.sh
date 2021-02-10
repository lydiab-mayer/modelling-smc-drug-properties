#!/bin/bash
#SBATCH --job-name=plotting
#SBATCH --account=smith
#SBATCH -o /scicore/home/smith/golmon00/MMC/TPP/JOB_OUT/plot.out
#SBATCH --mem=1GB
#SBATCH --qos=30min
#SBATCH --cpus-per-task=1
##################################
# Main script for plotting
#
# created 12.12.2019
# monica.golumbeanu@unibas.ch
#################################

SIM_FOLDER=$1
TEST_FOLDER=$2
PLOT_TITLE=$3

cd "plotting/"
# Plot the GP performance
echo "Plotting GP performance"
Rscript run_plot_performance.R $SIM_FOLDER $TEST_FOLDER 4 $PLOT_TITLE
Rscript run_plot_performance.R $SIM_FOLDER $TEST_FOLDER 6 $PLOT_TITLE

# Plot the sensitivity analysis
echo "Plotting sensitivity"
Rscript run_plot_sensitivity.R $SIM_FOLDER 4 $PLOT_TITLE
Rscript run_plot_sensitivity.R $SIM_FOLDER 6 $PLOT_TITLE

# Plot the optimization results



