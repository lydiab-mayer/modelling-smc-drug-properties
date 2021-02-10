
Halflife <- c(30,45,60,75,90,105,120,135,150)
Efficacy <- seq(0.7,1,by=0.05)
opt_var_name <- c("Halflife","Efficacy")
opt_setups <- expand.grid(Halflife,Efficacy)
names(opt_setups) <- c("Halflife","Efficacy")

table_file <- "/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/optimisation/incred/opt_setup_file.txt"
write.table(opt_setups, table_file, sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)



Coverage_SMC <- c(0.4,0.5,0.6,0.7,0.8,0.9,1)
Halflife <- c(30,45,60,75,90,105,120,135,150)
Efficacy <- seq(0.7,1,by=0.05)
opt_setups <- expand.grid(Coverage_SMC,Halflife,Efficacy)
names(opt_setups) <- c("Coverage_SMC_Res","Halflife","Efficacy")

table_file <- "/scicore/home/smith/burlyd00/smc_lai/E3E4_comp/gp_5/optimisation/non_inf/opt_setup_file.txt"
write.table(opt_setups, table_file, sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)


library(tgp)
rm(list = ls())

set.seed(42)

EIR <- c(0,0)
Coverage = c(0.4, 1)
Halflife =c(30,150)
Efficacy= c(0.7,1)
Coverage_SMC_Res = c(0.4,1)
param_ranges = rbind(EIR,Coverage,Coverage_SMC_Res,Halflife,Efficacy)

Xcand = as.data.frame(lhs(1000, param_ranges))


colnames(Xcand) = c("EIR","Coverage","Coverage_SMC","Halflife","Efficacy")

table_file <- "/scicore/home/penny/burlyd00/smc_lai/E3E4_comp/gp_5/optimisation/non_inf/opt_setup_file_lhs.txt"
write.table(Xcand, table_file, sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)

