#!/bin/bash
PARAM_TABLE_FILE=/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/param_tab_4.txt

SCAFFOLD_FILE=/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/scaffold.xml
BASE_FOLDER=/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/base_4/
SCENARIOS_FOLDER=/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/scenarios_4/
OM_FOLDER=/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/om/
NUM=$(wc -l < $PARAM_TABLE_FILE)
 echo "Creating $NUM-1 scenarios ..." 
#sbatch -W ../job_create_scenarios.sh $PARAM_TABLE_FILE $SCAFFOLD_FILE $BASE_FOLDER $SCENARIOS_FOLDER

 echo "Running OM simulations... " 
sbatch -W --array=1-$NUM ../run_OM.sh $SCENARIOS_FOLDER $OM_FOLDER

