# plot non-inferiority results
library(ggpubr)
library(dplyr)
library(ggplot2)
library(cowplot)   
library(patchwork) 

rm(list = ls())

mainDir  <- "~/smc_lai/analysis_workflow/analysis_scripts"
setwd(mainDir)
source('~/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R')
source('~/smc_lai/analysis_workflow/analysis_scripts/supp/E5_resources.R')

dir.create(file.path(paste0(mainDir,"/Outputs/E5_Compare/")), showWarnings = FALSE)


library(hetGP)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tgp)


SMC_HL <- c(31.295,20)
LAI_dec <- c("hill","exp","wei")
seasonaility <- c("Sen","Mali")
settings <- expand.grid(SMC_HL,LAI_dec,seasonaility)
names(settings) <- c("SMC_HL","LAI_dec","seasonality")



for (k in 1:nrow(settings)) {
  setting <- settings[k, ]
  train_GP_E5_noninf(setting)
}



toplot <-which(settings$LAI_dec=="hill" )


plotlist <- list()
for (k in 1:length(toplot) ) {
  
  setting <- settings[toplot[k], ]
  plotlist[[k]] <- pred_GP_E5_noninf(setting)
}


q1 <- ggarrange(plotlist = list(plotlist[[1]],
                               ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                               plotlist[[3]],
                               ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")),
               ncol=4,common.legend=FALSE,legend="bottom",heights=c(rep(1,4)), 
               widths=c(1,0.1,1,0.1), labels=c("a", "", "b",""),
               font.label = list(size = 40, face = "bold"),hjust=-0.1,vjust=1)

q12 <- ggarrange(plotlist = list(ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                q1  ),
                nrow=2,heights=c(c(0.05,1)))
q1an <-   annotate_figure(q12,
                           top = text_grob(  "No SP resistance" ,
                                             
                                             color = "black", face = "bold", size = 35),
                           
)
q2 <- ggarrange(plotlist = list(plotlist[[2]],
                               ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                               plotlist[[4]],
                               ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")),
               ncol=4,common.legend=FALSE,legend="bottom",heights=c(rep(1,4)), 
               widths=c(1,0.1,1,0.1), labels=c("c", "", "d",""),
               font.label = list(size = 40, face = "bold"),hjust=-0.1,vjust=1)

q22 <- ggarrange(plotlist = list(ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                 q2  ),
                 nrow=2,heights=c(c(0.05,1)))
q2an <-   annotate_figure(q22,
                          top = text_grob(  "Prevalent SP resistance" ,
                                            
                                            color = "black", face = "bold", size = 35),
                          
)

q3 <- ggarrange(plotlist = list(
  q1an,
  
                               q2an
                               ),
                nrow=2, 
                heights=c(1,1))

ggsave(
  filename =  paste0("./Outputs/E5_Compare/noninf_normSMC_test.jpg"),
  plot = last_plot(),
  
  width = 16,
  height = 16,
  dpi = 600)



for (k in 1:nrow(settings)) {
  setting <- settings[k, ]
  train_GP_E5_incred(setting)
}


toplot <-which(settings$SMC_HL>30)



plotlist <- list()
for (k in 1:length(toplot) ) {
  
  setting <- settings[toplot[k], ]
  plotlist[[k]] <- pred_GP_E5_incred(setting)
}


q <- ggarrange(plotlist = list(
  ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
   plotlist[[2]][[1]],ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),plotlist[[5]][[1]],ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
  ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
  plotlist[[3]][[1]],ggparagraph(text=" ", face = "italic", size =0.5, color = "black"),plotlist[[6]][[1]],ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
  ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
   plotlist[[1]][[1]],ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),plotlist[[4]][[1]],ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
),
heights=c(0.1,1,0.1,1,0.1,1), 
widths=c(1,0.1,1,0.1
  
), 
labels=c(
  "", "", "","",
  
  "a", "", "b","",
  "", "", "","",
  
  "c","","d","",
  "", "", "","",
  
  "e","","f",""),
nrow=6,ncol=4,
font.label = list(size = 40, face = "bold"),hjust=-0.1,vjust=-0.2)

q1 <- ggarrange(plotlist = list(ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                plotlist[[1]][[2]],ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")) ,
                heights=c(0.5,1,0.5), 
                widths=c(1,1,1
                ),
                nrow=3,ncol=1)

q2 <- ggarrange(plotlist = list(q,ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                q1,ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")) ,
                heights=c(1,1,1,1), 
                widths=c(1.2,0.01,0.4,0.1
                ),
                nrow=1,ncol=4)




ggsave(
  filename =  paste0("./Outputs/E5_Compare/inc_red_normSMCtestlegend.jpg"),
  plot = last_plot(),
  
  width = 20,
  height = 21,
  dpi = 600,
  
)
