Basic workflow for creating scenario files and running OpenMalaria simulations

created by Monica Golumbeanu
monica.golumbeanu@unibas.ch
28.01.2019
__________________________________


DESCRIPTION:
The present workflow can be used to create scenario .xml files and run OpenMalaria simulations for each scenario.
The main script of the workflow is OM_base_workflow.sh and calls sequentially two array job scripts:
submit_create_scenario.sh which creates all the scenarios
run_OM.sh which runs OpenMalaria simulations for all the scenarios

Before running the workflow, make sure you modify the headers of the array job scripts such that the
resource requirements are adequate and the job diagnostic files are saved in the correct location by 
replacing the corresponding options in the job specification header:

#!/bin/bash
#SBATCH --job-name=run_RCD
#SBATCH --account=smith
#SBATCH -o /scicore/home/smith/golmon00/JOB_OUT/%A_%a.out
#SBATCH --mem=1G
#SBATCH --qos=30min
#SBATCH --cpus-per-task=1

REQUIRED RESOURCES:
OpenMalaria resource files needed for running the workflow are located in the folder OM_schema38/

You can change the OpenMalaria version by modifying the load module command in the script run_OM.sh
In addition, you will also need to copy the corresponding files as provided in OM_schema38/ for the new version
You can specify a new location for the resource folder in the script run_OM.sh

INPUT:
To run the workflow, first make sure you have an input folder containing the following two files:
scaffold.xml: xml file containing marked parameters (surrounded by @) which will be varied across simulations. 
param_tab.txt: text file with a parameter table where each row corresponds to one scenario to be generated.

The folder example/ provided with the workflow is and example of input folder containing a dummy parameter 
table (generated with the R script generate_param_table.R) and a scaffold.xml file. 
IMPORTANT: If you want to run the workflow on the example, please copy this folder to your 
local home directory. The workflow will create additional folders with outputs in this folder.

COMMAND TO RUN THE WORKFLOW:
bash OM_basic_workflow.sh input_folder
where input_folder should be replaced with a valid folder path.


OUTPUTS:
The workflow creates the following sub-folders located in the provided input_folder:
base/ contains for each scenario replacement patterns for the parameters 
scenarios/ contains the final scenario.xml files to be used for openMalaria simulations
om/ contains the results of the OpenMalaria simulations




