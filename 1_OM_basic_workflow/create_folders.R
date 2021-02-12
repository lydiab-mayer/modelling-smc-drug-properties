########################################
# script create_folders.R
#
# creates folders for storing simulation outputs
#
########################################

create_folders <- function(exp){
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  GROUP = "/scicore/home/penny/GROUP/M3TPP/"
  
  dir.create(paste0(GROUP,exp))
  dir.create(paste0("./Experiments/",exp,"/JOB_OUT"))
  dir.create(paste0("./Experiments/",exp,"/OM_JOBS"))
  
}