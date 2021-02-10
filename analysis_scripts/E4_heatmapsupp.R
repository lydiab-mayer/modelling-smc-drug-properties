


rm(list = ls())
library(plyr)
library(ggplot2)
library(coda)
library(ggpubr)
library(pracma)
library(gridExtra)
library(reshape2)
library(hetGP)
library(data.table)
library(Rsolnp)


load(paste0("./Outputs/E4_Fig4/heatmapsupp.RData"))
LAI_dec <- c("exp","wei","hill")
Access <- c(0.1,0.5)

plots <- expand.grid(LAI_dec,Access)
colnames(plots) <- c("LAI_dec","Access")

plotlist <- list()
for (k in 1:nrow(plots)) {
  toplot <- plots[k, ]
  
  
plot1 <- subset(alldf, Access==toplot[,"Access"] & LAI_dec ==toplot[,"LAI_dec"])


if(unique(plot1$Access)==0.5){ 
  labels <-     c(0.4,1,1.6,2.2,2.9)
  title <- "high access to healthcare"
}else{
  labels <-  c(0.45,  0.72,1.3,  2.3, 3.2) 
  title <- "low access to healthcare"
  plot1$optimal_lai_coverage <- ifelse(plot1$EIR==150,-1,plot1$optimal_lai_coverage)
}
plot1$EIR <- factor(plot1$EIR, levels = c(paste(unique(plot1$EIR))),
                    labels = labels)





p <- ggplot(data = plot1) +
  geom_tile(aes(x = Efficacy*100, y = Halflife, fill = optimal_lai_coverage*100)) +
  
  facet_grid(vars(regimen),vars(EIR))+
  scale_fill_gradientn(limits = c(40,100),
                       colours=c("#6d1c68",  "#f7921e"),
                       na.value = "grey") +
  labs(y = "Half-life of protective efficacy [days]", x = "Initial protective efficacy [%]",
       fill= "minimal \nLAI coverage [%] \nto ensure non-inferiority") +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x =element_text(size=15,colour="black",vjust=0.7),
        axis.text.y =element_text(size=15,colour="black"),
        axis.title=element_text(size=15,colour="black", face="bold"))+
  theme(strip.background = element_rect(colour="white", fill="white",
                                        size=1.5, linetype="solid"))+
  theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+ 
  theme(legend.text = element_text( size=15),
        legend.title= element_text(size=15, face="bold"))+   theme(legend.position="bottom")+
  
  theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+ guides(fill = guide_colourbar(barwidth = 13, barheight = 1.5))

p2 <- annotate_figure(p,   right = text_grob(unique(plot1$decay) , rot = 270,size=30,face="bold"),
                         top = text_grob(  expression(bold(paste("cpppy"["0.25-5y"]))) ,size=20,face="bold"))
plotlist[[k]]   <- p2           
}

pA <- ggarrange(plotlist=list(plotlist[[1]]+ theme(plot.background = element_rect(color = "black",size = 2))+
                                theme(plot.margin=unit(c(1,1,2,3),"pt")),
                              plotlist[[2]]+ theme(plot.background = element_rect(color = "black",size = 2))+
                                theme(plot.margin=unit(c(1,1,2,3),"pt")),
                              plotlist[[3]]+ theme(plot.background = element_rect(color = "black",size =2))+
                  theme(plot.margin=unit(c(1,1,2,3),"pt"))),
ncol = 1, nrow = 3,common.legend = TRUE, labels = c("a","b","c"),
font.label = list(size = 40, face = "bold"),legend="bottom")

ggsave(
  filename =  "./Outputs/E4_Fig4/pA_supp.pdf",
  plot = last_plot(),
  width = 12,
  height = 18,
  dpi = 2400
)  


pB <- ggarrange(plotlist=list(plotlist[[4]]+ theme(plot.background = element_rect(color = "black",size = 2))+
                                theme(plot.margin=unit(c(1,1,2,3),"pt")),
                              plotlist[[5]]+ theme(plot.background = element_rect(color = "black",size = 2))+
                                theme(plot.margin=unit(c(1,1,2,3),"pt")),
                              plotlist[[6]]+ theme(plot.background = element_rect(color = "black",size =2))+
                                theme(plot.margin=unit(c(1,1,2,3),"pt"))),
                ncol = 1, nrow = 3,common.legend = TRUE, labels = c("a","b","c"),
                font.label = list(size = 40, face = "bold"),legend="bottom")

ggsave(
  filename =  "./Outputs/E4_Fig4/pB_supp.pdf",
  plot = last_plot(),
  width = 12,
  height = 18,
  dpi = 2400
)  
