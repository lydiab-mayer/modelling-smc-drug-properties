#!/bin/bash
#


sbatch -W --array=1-1000 launch_optimise_coverage.sh 
