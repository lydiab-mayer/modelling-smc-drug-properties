SIM_FOLDER=$1
PREDICTED=$2

PARAM_RANGES_FILE=$SIM_FOLDER"param_ranges.RData"
GP_FOLDER=$SIM_FOLDER"gp/trained/"$PREDICTED"/"
OPT_DEST_DIR=$SIM_FOLDER"gp/optimisation/"$PREDICTED"/"
OPT_SETUP_FILE=$SIM_FOLDER"gp/optimisation/"$PREDICTED"/opt_setup.txt"


# Submit GP optimization analysis array job
gp_files=(${GP_FOLDER}*.RData)
NUM=${#gp_files[@]}

OPT_SETUP=$(wc -l < $OPT_SETUP_FILE)

for row_n in `seq 1 1 $OPT_SETUP`
do
    sbatch --array=1-$NUM job_optimise_parameter.sh $GP_FOLDER $PARAM_RANGES_FILE $OPT_DEST_DIR $OPT_SETUP_FILE $row_n
done
