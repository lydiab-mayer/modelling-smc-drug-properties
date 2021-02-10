source('~/smc_lai/analysis_workflow/analysis_scripts/supp/import_functions.R')

settings <- c("Sen","Mali")

plotlist <- list()
for (k in 1:length(settings )){

  setting <- settings[k]




# seasons = read.table("~/smc_lai/analysis_workflow/resource_files/seasonality.txt", sep="\t", header = TRUE)
# plot_Mali = as.data.frame(t(seasons[which(seasons$Seasonality == setting), 2:ncol(seasons)]))
# names(plot_Mali) <- "EIR"
# plot_Mali$time <- c(2, 8 ,14, 20, 26, 32, 38 ,44, 50, 56 ,62, 68)


 i=20
  param_table <- read.table(paste0("/scicore/home/penny/GROUP/smc_lai/E0_Cont/postprocessing_5/split/E0_Cont_",setting,"_4.9167_exp_0.1.txt"), 
                      header = T, as.is = TRUE, stringsAsFactors = FALSE)
  

  om_results_folder = paste0("/scicore/home/penny/GROUP/smc_lai/E0_Cont/om/")
  
  res_df <- list()

    #  param_table <- param_table[order(param_table$Scenario_Name), ]
    OM_result_file = paste(om_results_folder, param_table[i,]$Scenario_Name, "_",
                           +                                param_table[i,]$SEED, "_cts.txt", sep="")
    OM_result = read.table(OM_result_file, sep="\t",header=TRUE)
    
    res <- subset(OM_result, timestep %in% seq(3650-73,3650))
     res$plottime <- res$timestep-min(res$timestep)+1


param_table[i,c("Efficacy","Halflife","Coverage")]

prev <- import_prev("E0_Cont",i,param_table_file=paste0("/scicore/home/penny/GROUP/smc_lai/E0_Cont/postprocessing_5/split/E0_Cont_",setting,"_4.9167_exp_0.1.txt"))
prev <- subset(prev, year==10 & age_group %in% c(3))
prev$plottime <- prev$time-min(prev$time)+1
prev$norm <- prev$prev/sum(prev$prev)
#prev$time <- prev$time-min(prev$time)
to_plot <- subset(res, plottime %in% seq(1,73))
to_plot$norm <-  to_plot$simulated.EIR/sum(to_plot$simulated.EIR,)


if (setting=="Mali"){
  timingInt=c(47,53,59,65)
  int_timing = data.frame(timingInt=timingInt,
                          EIR= to_plot[timingInt,c("norm")])
  }else {
    timingInt=c(53,59,65)
    int_timing = data.frame(timingInt=timingInt,
                            EIR= to_plot[timingInt,c("norm")])  }


arrow.length <- 0.02
touchoff.distance <- 0.02 # distance between data and start of arrow
arrowhead.size <- 6 # in millimeters

if(setting=="Sen") {
  seas <- "short season"
} else {
  seas <- "long season"
}

p <- ggplot()+
geom_line(data=to_plot,aes(x=plottime,y=norm,linetype="EIR"),size=1)+
  geom_line(data=prev,aes(plottime,norm,linetype="prev"),size=1)+
  labs(x=" ",y="Proportion")+
  geom_segment(data = int_timing, aes(x = timingInt, y = EIR + touchoff.distance,
                                xend = timingInt, yend = EIR + touchoff.distance + arrow.length),
               arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "first", type = "closed"), 
               colour="#42733C", size = 2) +
  geom_segment(data = int_timing[1,], aes(x = timingInt, y = EIR + touchoff.distance+0.03,
                                      xend = timingInt, yend = EIR + touchoff.distance + arrow.length+0.03),
               arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "first", type = "closed"), 
               colour="#517594", size = 2) +
scale_linetype_discrete(name = " ",
                        labels = c(expression(paste("EIR")),
                        expression(paste(italic("Pf"), "PR"["2-10y"]))))+
   scale_x_continuous(
    breaks = seq(3,72,6),
    label = c("J","F","M","A","M","J","J","A","S","O","N","D")
  )+
  
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=15),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20,colour="black"),
        axis.text.x=element_text(margin = margin(t = 5)),
        axis.title=element_text(size=20,face="bold"), 
        legend.title=element_text(size=20,face="bold"))+ theme(legend.key=element_blank()) +
  theme(legend.background=element_blank()) + theme(legend.text = element_text(size=25,face="bold"))+
  theme(legend.key.width = unit(1.5,"cm"),legend.position = "right")+
  ggtitle(paste0(seas))+ylim(0,0.15)+
  theme(plot.title = element_text(size=22,face="bold"))

plotlist[[k]] <- p

}

p <- ggarrange(plotlist = plotlist, common.legend = TRUE, legend="bottom")

ggsave(
  filename= paste0("./Outputs/E0_Cont/EIRprevint.jpg"),
  plot = last_plot(),
  width = 10,
  height = 5,
  dpi = 1200)




sum(to_plot$simulated.EIR)
sum(to_plot$input.EIR)
which(to_plot$input.EIR==max(to_plot$input.EIR))
which(to_plot$simulated.EIR==max(to_plot$simulated.EIR))

seasons = read.table("~/smc_lai/analysis_workflow/resource_files/seasonality.txt", sep="\t", header = TRUE)
plot_Mali = as.data.frame(t(seasons[which(seasons$Seasonality == "Mali"), 2:ncol(seasons)]))
names(plot_Mali) <- "EIR"
plot_Mali$time <- c(2, 8 ,14, 20, 26, 32, 38 ,44, 50, 56 ,62, 68)+6
plot_Mali$EIRsc <-   plot_Mali$EIR/max(plot_Mali$EIR)

res <- import_EIRs("E5_2_CT_LAI",21)
prev <- import_prev("E5_2_CT_LAI",21)
prev <- subset(prev, year==1 & age_group %in% c(1,2,3))
to_plot <- subset(res, timestep %in% seq(2,2+72))

ggplot()+geom_point(data=to_plot,aes(x=timestep,y=input.EIR,colour="input"))+
  geom_point(data=to_plot,aes(x=timestep,y=simulated.EIR,colour="simulated"))+
  geom_point(data=plot_Mali,aes(time,EIRsc,colour="points"))+
  geom_point(data=prev,aes(time,prev,colour="prev"))+
  labs(x="time",y="EIR")+
  ggtitle(paste0("EIR= ",unique(round(to_plot$EIR) )) )
sum(to_plot$simulated.EIR)
sum(to_plot$input.EIR)
which(to_plot$input.EIR==max(to_plot$input.EIR))
which(to_plot$simulated.EIR==max(to_plot$simulated.EIR))





res <- import_EIRs("E4_SMCLAI",5571)
to_plot <- subset(res, timestep>3580)

ggplot(data=to_plot)+geom_point(aes(x=timestep,y=input.EIR,colour="input"))+
  geom_point(aes(x=timestep,y=simulated.EIR,colour="simulated"))+labs(x="time",y="EIR")+
  ggtitle(paste0("EIR= ",unique(round(to_plot$EIR) )) )

sum(to_plot$simulated.EIR)
sum(to_plot$input.EIR)
which(to_plot$input.EIR==max(to_plot$input.EIR))
which(to_plot$simulated.EIR==max(to_plot$simulated.EIR))

res <- import_EIRs("E4_SMCLAI",1)
to_plot <- subset(res, timestep>3580)

ggplot(data=to_plot)+geom_point(aes(x=timestep,y=input.EIR,colour="input"))+
  geom_point(aes(x=timestep,y=simulated.EIR,colour="simulated"))+labs(x="time",y="EIR")+
  ggtitle(paste0("EIR= ",unique(round(to_plot$EIR) )) )
sum(to_plot$simulated.EIR)
sum(to_plot$input.EIR)
which(to_plot$input.EIR==max(to_plot$input.EIR))
which(to_plot$simulated.EIR==max(to_plot$simulated.EIR))

res <- import_EIRs("E4_SMCLAI",51)
to_plot <- subset(res, timestep>3580)

ggplot(data=to_plot)+geom_point(aes(x=timestep,y=input.EIR,colour="input"))+
  geom_point(aes(x=timestep,y=simulated.EIR,colour="simulated"))+labs(x="time",y="EIR")+
  ggtitle(paste0("EIR= ",unique(round(to_plot$EIR) )) )
sum(to_plot$simulated.EIR)
sum(to_plot$input.EIR)
which(to_plot$input.EIR==max(to_plot$input.EIR))
which(to_plot$simulated.EIR==max(to_plot$simulated.EIR))
