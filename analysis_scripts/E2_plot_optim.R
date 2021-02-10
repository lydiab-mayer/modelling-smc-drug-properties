library(ggplot2)
setwd( "/scicore/home/smith/burlyd00/smc_lai/analysis_workflow/analysis_scripts/")

results_folder = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/optimisation/incred/"
# Load optim results 
table_file <- "/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/optimisation/incred/opt_setup_file.txt"
opt_var_name <- "Coverage"

opt_settings <- read.table( table_file,  header = TRUE)

EIRValuesRange = c(1,5,10,20,50,100,150,200)
inc_ranges = seq(0.1,1,0.1)

Access = c(0.1,0.5)
Seasonality <- c("Sen","Mali")
LAI_dec <- c("hill","exp")
IntAge <- c(4.9167)
scenarios <- expand.grid(Access,Seasonality,LAI_dec,IntAge)
names(scenarios) <- c("Access","Seasonality","LAI_dec","IntAge")


for (j in 1:nrow(scenarios)) {
scenario <- scenarios[j,]

opt_res <- list()
for (k in 1:nrow(opt_settings)){
  
  opt_file = paste(results_folder, "seeds_E2_LAI_",scenario[1,"Seasonality"],"_",scenario[1,"IntAge"],"_",
                   scenario[1,"LAI_dec"],"_",scenario[1,"Access"],"_incred_",opt_var_name,"_",k, ".RData", sep="")
  
  load(opt_file)
  
 opt_res[[k]] <-  data.frame(cbind(opt_obj$table ,opt_obj$opt_setup))
 }
  
 res <- cbind( bind_rows(opt_res))
p1 <- ggplot(data = res) + theme(legend.position="bottom")+
  geom_tile(aes(x = Halflife, y = Efficacy, fill = opt_param)) +
  facet_grid( incred~ EIR ) + #, labeller = labeller(EIR_band = EIR.labs, SMC_coverage_band = SMC.labs)) +
  scale_fill_gradientn(limits = c(0.4,1),
                       colours=c("navyblue",  "darkorange1")) +
  labs( fill= "LAI coverage" )+theme(axis.title = element_text(size = 15))

annotate_figure(p1,
                top = text_grob("EIR", face = "bold", size = 14),
               
                
                right = text_grob("Incidence reduction", rot = 270),
                fig.lab =  paste0("Optimal LAI coverage: Seas:",scenario[1,"Seasonality"]," IntAge:",scenario[1,"IntAge"]," LAI decay:",
                                                                                     scenario[1,"LAI_dec"]," Access:",scenario[1,"Access"]), fig.lab.face = "bold", fig.lab.size = 20)

save_id <- paste(scenario[1,"Seasonality"],"_",scenario[1,"IntAge"],"_",
                 scenario[1,"LAI_dec"],"_",scenario[1,"Access"], sep="")
save_plot(
  filename =  paste0("./Outputs/E2_Optim/Optim_Coverage_",save_id,".jpg"),
  fig = last_plot(),
  width = 25,
  height = 20,
  dpi = 600
)  
}