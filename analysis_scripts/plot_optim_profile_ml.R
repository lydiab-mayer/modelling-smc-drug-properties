library(ggplot2)
setwd( "/scicore/home/smith/laagmi01/smc_lai/analysis_workflow/analysis_scripts/")


# Load optim results 
load('outfiles/test_ml.RData')

scenarios <- opt_res$scenarios

n_EIR_bands <- 3
EIR_steps <- (max(scenarios$EIR) - min(scenarios$EIR))/n_EIR_bands

scenarios$EIR_band <- 0
for (k in 1:n_EIR_bands){
  lower <- min(scenarios$EIR) + (k-1)*EIR_steps
  upper <- min(scenarios$EIR) + k*EIR_steps
  index <- which(scenarios$EIR >= lower & scenarios$EIR <= upper)
  scenarios[index, ]$EIR_band <- k
}

n_SMC_coverage_bands <- 3
SMC_coverage_steps <- (max(scenarios$Coverage_SMC) - min(scenarios$Coverage_SMC))/n_SMC_coverage_bands

scenarios$SMC_coverage_band <- 0
for (k in 1:n_SMC_coverage_bands){
  lower <- min(scenarios$Coverage_SMC) + (k-1)*SMC_coverage_steps
  upper <- min(scenarios$Coverage_SMC) + k*SMC_coverage_steps
  index <- which(scenarios$Coverage_SMC >= lower & scenarios$Coverage_SMC <= upper)
  scenarios[index, ]$SMC_coverage_band <- k
}

EIR.labs <- c("low EIR (1 - 9)", "medium EIR (9 - 17) ","high EIR (17 - 25)")
names(EIR.labs) <- c(1,2,3)

SMC.labs <- c("low SMC coverage (0.4 - 0.6)", "medium SMC coverage (0.6 - 0.8)","high SMC coverage (0.8 - 1.0)")
names(SMC.labs) <- c(1,2,3)


p <- ggplot(data = scenarios) +
  geom_tile(aes(x = Halflife, y = Efficacy, fill = optimal_lai_coverage)) +
  facet_grid( EIR_band ~ SMC_coverage_band, labeller = labeller(EIR_band = EIR.labs, SMC_coverage_band = SMC.labs)) +
  scale_fill_gradientn(limits = c(0.5,1),
                       colours=c("navyblue",  "darkorange1")) +
  labs(title = "optimal LAI coverage", subtitle = "optimal = lowest possible coverage that achieves 0.8 probability of non-inferiority", fill= "LAI coverage")


ggsave('outfiles/plot_optimisation.pdf', p, width = 15, height = 15)
