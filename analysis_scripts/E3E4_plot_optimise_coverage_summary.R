

rm(list = ls())
setwd("/scicore/home/penny/burlyd00/smc_lai/analysis_workflow/analysis_scripts")
library(nloptr)
library(hetGP)
library(Rsolnp)
library(metaheuristicOpt)
library(tools)
library(rapportools)
library(stringr)
library(ggplot2)
library(ggpubr)

ranges_file = "/scicore/home/smith/GROUP/smc_lai/E4_SMCLAI/param_ranges.RData"
load(ranges_file)
rownames(param_ranges)[3] <- "Coverage_SMC_Res"

results_folder = "/scicore/home/penny/burlyd00/smc_lai/E3E4_comp/gp_5/optimisation/non_inf/"
# Load optim results 
table_file <- "/scicore/home/penny/burlyd00/smc_lai/E3E4_comp/gp_5/optimisation/non_inf/opt_setup_file_lhs.txt"
opt_var_name <- "Coverage"

opt_settings <- read.table( table_file,  header = TRUE)




Access <- c(0.1,0.5)
LAI_dec <- c("hill","exp")
seasonality <- c("Mali","Sen")
IntAge <- c(4.9167)
settings <- expand.grid(Access,LAI_dec,seasonality,IntAge)
names(settings) <- c("Access","LAI_dec","seasonality","IntAge")

i=3
setting <- settings[i,]


test <- read.table( file=paste0("/scicore/home/penny/burlyd00/smc_lai/E3E4_comp/gp_5/optimisation/processed_ML/scenarios_Mali_4.9167_hill_0.1.txt"),  header = TRUE)
test1 <- test[test$Halflife==110 & test$Efficacy==0.7 , ]

test2 <- test[which(test$Coverage_SMC_Res>test$optimal_lai_coverage & test$optimal_lai_coverage>0) ,]





dfm = plyr::rename(res,c("Coverage_SMC_Res"="x","Halflife"="y","Efficacy"="z"))
crd = coord_tern()
rev = xy2tlr(dfm,crd)

plot <- res[which(res$EIR==11),]
plot <- plot[c(1,100,200),]
ggtern(test2,aes(x=Efficacy,y=Coverage_SMC_Res,z=Halflife,color=EIR)) +
   theme_rgbw() + 
  geom_point()+                          #define data sources
  scale_color_viridis_c()+
  theme_showarrows() +theme_clockwise() + labs( x       = "Efficacy",
                                                y       = "SMC_Coverage",
                                                z       = "Halflife")+
annotate(geom  = 'text',
           x     = c(0,0,0),
           y     = c(0.5,1/3,0.0),
           z     = c(0.0,1/3,1.0),
           angle = c(0,30,60),
           vjust = c(1.5,0.5,-0.5),
           label = paste("Point",c("A","B","C")),
           color = c("green","red",'blue')) 



plot_opt_E34 <- function(setting) {
  
 
  opt_res <- list()
  for (k in 1:nrow(opt_settings)){
    
    opt_file = paste(results_folder, "seeds_E3E4_comp_",setting[1,"seasonality"],"_",setting[1,"IntAge"],"_",
                     setting[1,"LAI_dec"],"_",setting[1,"Access"],"_UL_",opt_var_name,"_",k, ".RData", sep="")
    
    load(opt_file)
    df <- opt_obj$table
df$Coverage_SMC_Res <- opt_obj$opt_setup["Coverage_SMC_Res"]
  df$Halflife <- opt_obj$opt_setup["Halflife"]
  df$Efficacy <- opt_obj$opt_setup["Efficacy"]
  

opt_res[[k]] <-  df
  }
  
  res <-  as.data.frame(do.call(rbind, opt_res))
 
  res <- subset(res, EIR==5)

  test <- res[!complete.cases(res),]
 res <-  test2
res$Halflife <- 100* (res$Halflife-min(res$Halflife))/(max(res$Halflife)-min(res$Halflife))
res$Efficacy <-  100* (res$Efficacy-min(res$Efficacy))/(max(res$Efficacy)-min(res$Efficacy))
res$Coverage_SMC_Res <-  100* (res$Coverage_SMC_Res-min(res$Coverage_SMC_Res))/(max(res$Coverage_SMC_Res)-min(res$Coverage_SMC_Res))
#df$opt_param <-  100* (df$opt_param-min(df$opt_param))/(max(df$opt_param)-min(df$opt_param))
res$opt_param <- ifelse(res$opt_param==0,1.2,res$opt_param)
res$opt_param <- 100* (res$opt_param-min(res$opt_param))/(max(res$opt_param)-min(res$opt_param))


  
res <- res[complete.cases(res),]

  
  
  # res$opt_param <- ifelse(res$opt_param==0,0.4,res$opt_param)
  
 n_points = round(nrow(res)*0.9)
 
 res$id <- seq(1,nrow(res))
 index_train = sample(nrow(res), n_points)
 index_test = setdiff(1:nrow(res), index_train)
 
 train_data <-res[res$id %in% index_train,]
 test_data <- res[res$id %in% index_test,]

 # train GP on upper limit of the confidence interval
 predicted= "opt_param"
 predictors <- c("Efficacy","Halflife","Coverage_SMC_Res")
 
 train_model_matern <- train_GP_matern(train_data[ ,c(predictors, predicted)])
 
 
 test_GP_plot(GP_model=train_model_matern, test=test_data[,c(predictors, predicted)],
              save=paste(results_folder,"Coverage", "_R2.jpg", sep=""))
 
 
 
 Efficacy <- seq(0,100,1)
 Halflife <- seq(0,100,1)
 Coverage_SMC_Res <-  seq(0,100,1)
 
 Xcand <- expand.grid(Halflife,Efficacy,Coverage_SMC_Res)
 names(Xcand) <- c("Halflife","Efficacy","Coverage_SMC_Res")
 
 
 prediction_for_plot <- predict(x = as.matrix(Xcand), 
                                object = train_model_matern)
 upper_limit_CI   <- data.frame(cbind(Xcand,pred=prediction_for_plot$mean ))
 
 upper_limit_CI$pred <- ifelse( upper_limit_CI$pred>=83,NA, upper_limit_CI$pred)
 
 breaks = seq(0,100,10)
 
 ggtern(upper_limit_CI,aes(Halflife,Efficacy,Coverage_SMC_Res,colour=pred)) +
   geom_point()+                          #define data sources
   scale_color_viridis_c()+
   theme_showarrows() +theme_clockwise() + labs( x       = "Halflife",
                                                 y       = "Efficacy",
                                                 z       = "SMC_Coverage")
 
 

 
 ggtern(res,aes(Halflife,Efficacy,Coverage_SMC_Res,colour=EIR)) +
   geom_point()+                          #define data sources
   scale_color_viridis_c()+
   theme_showarrows() +theme_clockwise() + labs( x       = "Halflife",
                                                 y       = "Efficacy",
                                                 z       = "SMC_Coverage")
 
 
 
 
 
df <- res[res$EIR==11, ]  c("Efficacy","Halflife","Coverage_SMC_Res","optimal_lai_coverage")]
names(df) <- c("x","y","z","LAI_cover")
df$x <- 100* (df$x-min(df$x))/(max(df$x)-min(df$x))
df$y <-  100* (df$y-min(df$y))/(max(df$y)-min(df$y))
df$z <-  100* (df$z-min(df$z))/(max(df$z)-min(df$z))
df$LAI_cover <- ifelse(df$LAI_cover<0,NA,df$LAI_cover)
df$LAI_cover <-  100* (df$LAI_cover-min(df$LAI_cover))/(max(df$LAI_cover)-min(df$LAI_cover))

breaks = seq(0,100,10)
  
 ggtern(df,aes(x,y,z,colour=LAI_cover)) +
   geom_point()+                          #define data sources
   scale_color_viridis_c()+
   theme_showarrows() +theme_clockwise() + labs( x       = "Efficacy",
         y       = "Halflife",
         z       = "SMC_Coverage")+ 
   geom_interpolate_tern(aes(value=LAI_cover,color=..level..)) 

  ggtern(df,aes(x,y,z)) + 
   geom_interpolate_tern(aes(value= LAI_cover,fill=..level..),
                         colour   = "white",
                         formula  = value~poly(x,y,
                                               degree=2,
                                               raw=TRUE),
                         method   = "lm",
                         binwidth = 25,
                         buffer   = 1.5,
                         n        = 200)
 
   geom_interpolate_tern(aes(value=level),
                         base = 'identity',method = glm,
                         formula = value ~ polym(x,y,degree = 2),
                         n = 150, breaks = breaks)
 
   scale_T_continuous(breaks=unique(df$y),labels=labFnc(df$y)) +
   scale_L_continuous(breaks=unique(df$x),labels=labFnc(df$x)) +
   scale_R_continuous(breaks=unique(df$z),labels=labFnc(df$y)) +
   theme_bw() +
   theme_nogrid_minor() + 
   geom_point() +
   labs(title="Example Use of Data Along Ternary Axes")
   
   
  #define a data geometery with an aesthetic 
   
plots <- list()
for (l in 1:nrow(settings)){
  print(l)
  
  seasonality <- settings[l,1]
  decay <- settings[l,2]
  acess <- settings[l,3]
  age <- settings[l,4]
  
  scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess, sep = "")
  
  scenarios <- read.table(paste('scenarios_processed/scenarios_', scen_name, '.txt', sep = ""), header = T)
  
  
  scenarios_plot <- scenarios
  index <- which(scenarios$EIR == 11 | scenarios$EIR == 51)
  scenarios_plot <- scenarios_plot[index,]
  index <- which(scenarios_plot$Coverage_SMC_Res == 0.4 | scenarios_plot$Coverage_SMC_Res == 0.8)
  scenarios_plot <- scenarios_plot[index,]
  

  EIR.labs <- c("EIR = 11", "EIR = 51")
  names(EIR.labs) <- c(11,51)
  
  SMC.labs <- c("SMC coverage = 0.4", "SMC coverage = 0.8")
  names(SMC.labs) <- c(0.4,0.8)
  
  
  plots[[l]] <- ggplot(data = scenarios_plot) +
    geom_tile(aes(x = Halflife, y = Efficacy, fill = optimal_lai_coverage)) +
    facet_grid( EIR ~ Coverage_SMC_Res, labeller = labeller(EIR = EIR.labs, Coverage_SMC_Res = SMC.labs)) + 
    scale_fill_gradientn(limits = c(0.4,1),
                         colours=c("green4",  "grey"),
                         na.value = "white") +
    labs(title = paste("Seasonality: ", seasonality, ", Decay: ", decay, sep =""), fill= "LAI coverage")
  
  
 
  
}

p1 <- plots[[1]]
p2 <- plots[[2]]
p3 <- plots[[3]]
p4 <- plots[[4]]

p <- ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2)



ggsave(paste('outfiles/summary_plot_coverage.pdf', sep = ""), p, width = 15, height = 10)

