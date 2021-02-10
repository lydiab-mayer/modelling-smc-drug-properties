
# plot non-inferiority results
library(ggpubr)
library(dplyr)
library(grid)


rm(list = ls())

mainDir  <- "~/smc_lai/analysis_workflow/analysis_scripts"
setwd(mainDir)
source('~/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R')


dir.create(file.path(paste0(mainDir,"/Outputs/E5_Compare/")), showWarnings = FALSE)


library(hetGP)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tgp)


LAI_dec <- c("hill")
seasonality <- c("Sen","Mali")
Access <- c(0.1,0.5)
settings <- expand.grid(Access,LAI_dec,seasonality)
names(settings) <- c("Access","LAI_dec","seasonality")

plotlist <- list()
for (k in 1:nrow(settings)){
setting <- settings[k, ]



if(unique(setting$Access)==0.5){ 
  EIRs <- c(5,9,20,47,150)
  labels <-     c(0.4,1,1.6,2.2,2.9)
  title <- "high access to healthcare"
}else{
  EIRs <- c(3,4,8,28,150)
  labels <-  c(0.45,  0.72,1.3,  2.3, 3.2) 
  title <- "low access to healthcare"
}

if(setting$seasonality=="Sen") {
  seas <- "short season"
} else {
  seas <- "long season"
}
datalist <- list()
predlist <- list()
Cov60list <-list()
for (j in 1:length(EIRs)) {
  EIR <- EIRs[j]
scen_name <- paste(setting[,"seasonality"],"_4.9167_",setting[,"LAI_dec"],"_",setting[,"Access"],"_",EIR, sep = "")

outfile_pppy_SMC <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_SMC_', scen_name, '.RData', sep = "")

load(outfile_pppy_SMC)


Coverage_SMC <-  seq(0.4,1,0.01)

prediction_for_plot <- predict(x = as.matrix(Coverage_SMC), 
                               object = GP_trained_SMC)

pred   <- data.frame(cbind(Coverage_SMC,pred=prediction_for_plot$mean ))
pred$EIR <- EIR
pred$EIR_group <- j
pred$init <- labels[j]
predlist[[j]] <- pred

Cov60 <- (pred$init-pred$pred)/pred$init
Cov60diff <- Cov60[61]-Cov60[23]
Cov60list[[j]] <- data.frame(cbind(unique(pred[,c("EIR","EIR_group","init")]),Cov60diff ))
filename_SMC <- paste('/scicore/home/penny/GROUP/smc_lai/E3_SMCSMCdisc/postprocessing_5/seeds_E3SMCSMCdisc_', scen_name, '.txt', sep = "")

data <- read.table(filename_SMC, header = T)


filename_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/postprocessing_5/seeds_E4SMCLAIdisc_', scen_name, '.txt', sep = "")
test <- read.table(filename_LAI, header = T)

data$EIR_group <- j
data$init <-labels[j]

datalist[[j]] <- data
}

allpred <- rbindlist(Cov60list)

alldf <- rbindlist(predlist)


brewer.pal( n=9, name = "Greens")
colors <- c(  "1"="#E5F5E0", "2"="#A1D99B" ,"3"="#74C476","4" ="#238B45","5"= "#00441B")

p <- ggplot()+ geom_line(data= alldf, aes(x=Coverage_SMC*100 , y = 100*(init-pred)/init, 
                                         color = as.character(EIR_group)),size=1.5)  +
  geom_vline(aes(xintercept = 62), colour = "grey50", linetype = "dashed",size=1.5) +
  geom_vline(aes(xintercept = 100), colour = "grey50", linetype = "dashed",size=1.5) +
  geom_point(data= alldf[which(alldf$Coverage_SMC %in% c(0.62,1)),], aes(x=Coverage_SMC*100 , y = 100*(init-pred)/init, 
                             fill = as.character(EIR_group)),shape=21,size=3) +
  scale_fill_manual(values=colors,
                     labels=round(allpred$Cov60diff*100, digits = 0),
                     breaks=c("1","2","3","4","5"),drop = FALSE) +
  scale_color_manual(values=colors,
                     labels=labels,
                     breaks=c("1","2","3","4","5"),drop = FALSE) +
  labs(fill= expression(paste("additional \nincidence reduction"["0.25-5y"] )) ,color =   expression(paste("initial cases per \nperson per year"["0.25-5y"])), 
       x ="Coverage SMC-SP+AQ [%]", y = expression(paste("incidence reduction"["0.25-5y"],"[%]")  ) )+ylim(0,100)+
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20,colour="black"),
        axis.text.x=element_text(margin = margin(t = 5)),
        axis.title=element_text(size=20,face="bold"), 
        legend.title=element_text(size=20,face="bold"))+ theme(legend.key=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.text = element_text(size=15,face="bold"))+
  theme(legend.key.width = unit(1.5,"cm"))+ theme(legend.key=element_blank(),  legend.text = element_text(size=20,face="bold"))+
  theme(plot.title = element_text(size = 22, face = "bold",hjust = 0.2, vjust=2.12))+
  ggtitle(paste0(seas,"; ", title))+
  geom_segment( aes(x = 62, y =85, xend = 100, yend = 85),
                colour = "black",size=1.5,linejoin = c('bevel'),
                arrow = arrow(length = unit(0.3, "inches")) )+ theme(plot.background = element_rect(color = "black",size = 2))+
  theme(plot.margin=unit(c(15,10,10,10),"pt"))+ guides(fill = guide_legend(order = 1),
                                                    colour = guide_legend(order = 2))



plotlist[[k]] <- p
}

p1 <- ggarrange(plotlist = list(plotlist[[1]],plotlist[[3]]) ,labels = c("a","b"),
                font.label = list(size = 40, face = "bold"),hjust=-0.5,vjust=2)
p2 <- ggarrange(plotlist = list(plotlist[[2]],plotlist[[4]]), labels = c("c","d"),
                font.label = list(size = 40, face = "bold"),hjust=-0.5,vjust=2)

p3 <- ggarrange(plotlist = list( 
p1,p2),nrow=2)
ggsave(
  filename= paste0("./SMCpppy.jpg"),
  plot = last_plot(),
  
  width = 14,
  height = 13,
  dpi = 600 )



