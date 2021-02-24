#!/bin/bash
#
# gp_file = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/gp/trained/prevred_int_y10/seeds_E0MAB_Mali_4.9167_exp_0.1_10_prevred_int_y10_cv.RData"
# ranges_file = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/param_ranges.RData"
# results_folder = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/optimisation/"

gp_file=$1
ranges_file=$2
results_folder=$3

sbatch job_optimise_coverage.sh $gp_file $ranges_file $results_folder
