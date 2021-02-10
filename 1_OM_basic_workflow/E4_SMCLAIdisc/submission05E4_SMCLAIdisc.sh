#!/bin/bash
PARAM_TABLE_FILE=/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/param_tab_5.txt

SCAFFOLD_FILE=/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/scaffold.xml
BASE_FOLDER=/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/base_5/
SCENARIOS_FOLDER=/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/scenarios_5/
OM_FOLDER=/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/om/
NUM=$(wc -l < $PARAM_TABLE_FILE)
 echo "Creating $NUM-1 scenarios ..." 
sbatch -W ../job_create_scenarios_ml.sh $PARAM_TABLE_FILE $SCAFFOLD_FILE $BASE_FOLDER $SCENARIOS_FOLDER

 echo "Running OM simulations... " 
sbatch -W --array=1-$NUM ../run_OM_ml.sh $SCENARIOS_FOLDER $OM_FOLDER

