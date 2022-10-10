SIM_FOLDER=$1
PREDICTED=$2
OPT_SETUP_FILE=$3
n_gridpoints=$4
SCALE=$5

PARAM_RANGES_FILE=$SIM_FOLDER"param_ranges.RData"
GP_FOLDER=$SIM_FOLDER"gp/trained/"$PREDICTED"/"
OPT_DEST_DIR=$SIM_FOLDER"gp/optimisation/"$PREDICTED"/"

# Submit GP optimization analysis array job
gp_files=(${GP_FOLDER}*.RData)
NUM=${#gp_files[@]}

sbatch --array=1-$NUM --time=1-00:00:00 job_optimise_parameter.sh $GP_FOLDER $PARAM_RANGES_FILE $OPT_DEST_DIR $OPT_SETUP_FILE $n_gridpoints $SCALE

