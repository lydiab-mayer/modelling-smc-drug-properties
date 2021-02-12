########################################
# script generate_experiment_param_table.R
#
# construct a table with parameter values for each scenario of an OM experiment
#
# INPUTS:
#   table_file - name of the file where the table should be saved
#   experiment - name of the experiment
#
# OUTPUTS:
#	- a csv table with all the parameters for OM experiment

# created 21.01.2018
# monica.golumbeanu@unibas.ch
########################################
gen_paramtable <- function(exp, param_ranges_cont,param_cat, noSamples, noSeeds,chunk_size) {

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

GROUP = "/scicore/home/penny/GROUP/M3TPP/"

ranges_file = paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/param_ranges.RData")
table_file = paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/param_tab.txt")
scaffold_file= paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_jobs/scaffold.xml")


# sample continuous parameter values 

Xcand = as.data.frame(lhs(noSamples, param_ranges_cont))
colnames(Xcand) = rownames(param_ranges_cont)

SEED = c(1:noSeeds)

# save continous parameter ranges for the experiment and copy to GROUP folder
save(param_ranges_cont, file = ranges_file)

file.copy(ranges_file, paste0("../../GROUP/M3TPP/",exp,"/param_ranges.RData"),overwrite=TRUE)



# Table with the parameter values
param_all <- param_cat
param_all[[length(param_cat)+1]] <- as.data.frame(Xcand)

param_tab_all = Reduce(merge ,param_all)



param_tab = cbind(scenarios_names= paste(exp, 1:nrow(param_tab_all), sep="_"),
                   param_tab_all)
 colnames(param_tab)[1] = "Scenario_Name"
 
# Add seed column at the beginning
 param_tab = merge(param_tab, as.data.frame(SEED))
 
 # Write table to specified destination file
file.remove(table_file)
 write.table(param_tab, table_file, sep = "\t", quote = FALSE, col.names = TRUE,
             row.names = FALSE)
 
 file.copy(table_file, paste0("../../GROUP/M3TPP/",exp,"/param_tab.txt"),overwrite=TRUE)
 
 # copy over scaffold file 
 file.remove(paste0(GROUP,"/",exp,"/scaffold.xml"))

 file.copy(scaffold_file, paste0("../../GROUP/M3TPP/",exp,"/scaffold.xml"),overwrite=TRUE)




chunk <- chunk_size
n <- nrow(param_tab)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
split_param_tab <- split(param_tab,r)

num_tab <- max(r)
for(j in 0:(num_tab-1)){ 
  
table_file = paste0(GROUP,"/",exp,"/param_tab_",j,".txt")

# Write table to specified destination file
file.remove(table_file)
write.table(split_param_tab[[j+1]], table_file, sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)

}

}