#on cluster: ls >> filenames_of_simulations_done.txt

rm(list = ls())
library(sjmisc)
experiment <- "E5_2_CT_LAI"
f <- read.table(paste0("/scicore/home/penny/burlyd00/smc_lai/analysis_workflow/1_OM_basic_workflow/",experiment, "_simulations_done.txt"),header = F, stringsAsFactors = F) 
a1 <- read.table(paste0("/scicore/home/penny/GROUP/smc_lai/",experiment, "/param_tab_0.txt"),header = T, stringsAsFactors = F)
a2 <- read.table(paste0("/scicore/home/penny/GROUP/smc_lai/",experiment, "/param_tab_1.txt"),header = T, stringsAsFactors = F)
a3 <- read.table(paste0("/scicore/home/penny/GROUP/smc_lai/",experiment, "/param_tab_2.txt"),header = T, stringsAsFactors = F)
a4 <- read.table(paste0("/scicore/home/penny/GROUP/smc_lai/",experiment, "/param_tab_3.txt"),header = T, stringsAsFactors = F)
a5 <- read.table(paste0("/scicore/home/penny/GROUP/smc_lai/",experiment, "/param_tab_4.txt"),header = T, stringsAsFactors = F)
a6 <- read.table(paste0("/scicore/home/penny/GROUP/smc_lai/",experiment, "/param_tab_5.txt"),header = T, stringsAsFactors = F)

a <- data.frame(rbind(a1,a2,a3,a4,a5,a6))


f_1 <- f[grepl("_out.txt", f$V1), "V1"]

sims <- paste0(ldply(strsplit(f_1,"_"))[,5], "_",ldply(strsplit(f_1,"_"))[,6])
params <- paste0(ldply(strsplit(a$Scenario_Name,"_"))[,5], "_",a$SEED)

simulated <- which(params %in% sims)

param_tab_new <- a[-simulated,]


write.table(param_tab_new,paste0("/scicore/home/penny/GROUP/smc_lai/",experiment, "/param_tab_new.txt"), sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)