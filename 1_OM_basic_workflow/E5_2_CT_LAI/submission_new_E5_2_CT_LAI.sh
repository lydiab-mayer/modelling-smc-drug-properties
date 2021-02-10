#!/bin/bash
PARAM_TABLE_FILE=/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/param_tab_new.txt

SCAFFOLD_FILE=/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/scaffold.xml
BASE_FOLDER=/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/base_new/
SCENARIOS_FOLDER=/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/scenarios_new/
OM_FOLDER=/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/om/
NUM=$(wc -l < $PARAM_TABLE_FILE)
 echo "Creating $NUM-1 scenarios ..." 
sbatch -W ../job_create_scenarios.sh $PARAM_TABLE_FILE $SCAFFOLD_FILE $BASE_FOLDER $SCENARIOS_FOLDER

rm -r $BASE_FOLDER
 echo "Running OM simulations... " 
sbatch -W --array=1-$NUM ../run_OM.sh $SCENARIOS_FOLDER $OM_FOLDER

