#################################
# Plot the seasonality pattern
#
#
# created 03.10.2019
# monica.golumbeanu@unibas.ch
################################

library(data.table)

seasons = read.table("~/smc_lai/analysis_workflow/resource_files/seasonality.txt", sep="\t", header = TRUE)
plot_Mali = as.data.frame(t(seasons[which(seasons$Seasonality == "Mali"), 2:ncol(seasons)]))
plot_Sen = as.data.frame(t(seasons[which(seasons$Seasonality == "Sen"), 2:ncol(seasons)]))


plot_df <- plot_Sen
colnames(plot_df) = "Seasonality"
plot_df$Month = c(1:12)
plot_df$Seasonality = plot_df$Seasonality/sum(plot_df$Seasonality)

res1 = plot_df[c("Aug"),]
res2 = plot_df[c("Sep"),]
res3 = plot_df[c("Oct"),]
res4 = plot_df[c("Nov"),]

arrow.length <- 0.08
touchoff.distance <- 0.03 # distance between data and start of arrow
arrowhead.size <- 2 # in millimeters

p = ggplot(plot_df, aes(x=Month, y=Seasonality)) + geom_line() + geom_point() +
    ylim(0, 0.9) + labs( x = "Month", y = "Proportion", title = "Modelled monthly transmission pattern") + 
    scale_x_continuous(breaks = c(1:12), labels = c(rownames(plot_df))) +
    theme(strip.background = element_rect(colour="white", fill="white")) +
    geom_segment(data = res4, aes(x = Month, y = Seasonality + touchoff.distance,
                                xend = Month, yend = Seasonality + touchoff.distance + arrow.length),
                                arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "first"), 
                 colour="#2ca25f", size = 1) +
  geom_segment(data = res2, aes(x = Month, y = Seasonality + touchoff.distance,
                                xend = Month, yend = Seasonality + touchoff.distance + arrow.length),
               arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "first"), 
               colour="#2ca25f", size = 1) +
  geom_segment(data = res3, aes(x = Month, y = Seasonality + touchoff.distance,
                                xend = Month, yend = Seasonality + touchoff.distance + arrow.length),
               arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "first"), 
               colour="#2ca25f", size = 1) +
    geom_segment(data = res1, aes(x = Month, y = Seasonality + touchoff.distance+0.1,
                                xend = Month, yend = Seasonality + touchoff.distance +0.1+ arrow.length),
                                arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "first"),
                 colour="#2b8cbe", size = 1)+
  
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=13),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=13,colour="black"),
        axis.text.x=element_text(margin = margin(t = 5)),
        axis.title=element_text(size=13,face="bold"),  
        legend.title=element_text(size=13))+
  theme(legend.background=element_blank()) + theme(legend.key=element_blank())

plot(p)


# plot for presentation/paper 

plot_df2 <- plot_df[7:12,]

colors <- c("SMC" = "#2ca25f", "LAI"="#2b8cbe")

ggplot(plot_df2, aes(x=Month, y=Seasonality)) + geom_line() + geom_point() +
  theme_bw(base_size=13) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(0, 0.75) + labs( x = "Month", y = "% yearly EIR",color = "Intervention") + 
  scale_color_manual(values = colors)+
  
  scale_x_continuous(breaks = c(1:12), labels = c(rownames(plot_df))) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  geom_segment(data = res4, aes(x = Month, y = Seasonality + touchoff.distance,
                                xend = Month, yend = Seasonality + touchoff.distance + arrow.length, color="SMC"),
               arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "first"), size = 2,
               show.legend = TRUE) +
  geom_segment(data = res2, aes(x = Month, y = Seasonality + touchoff.distance,
                                xend = Month, yend = Seasonality + touchoff.distance + arrow.length, color="SMC"),
               arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "first"), 
                size = 2) +
  geom_segment(data = res3, aes(x = Month, y = Seasonality + touchoff.distance,
                                xend = Month, yend = Seasonality + touchoff.distance + arrow.length, color="SMC"),
               arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "first"), 
                size = 2) +
  geom_segment(data = res1, aes(x = Month, y = Seasonality + touchoff.distance+0.1,
                                xend = Month, yend = Seasonality + touchoff.distance +0.1+ arrow.length, color="LAI"),
               arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "first"),
               size = 2)+
   theme(legend.position = c(0.8, 0.8))+
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=13),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=13,colour="black"),
        axis.text.x=element_text(margin = margin(t = 5)),
        axis.title=element_text(size=13,face="bold"),  
        legend.title=element_text(size=13))+
  theme(legend.background=element_blank()) + theme(legend.key=element_blank())


save_plot(
  filename =  paste0( "./Outputs/Transmission_Mali.jpg"),
  fig = last_plot(),
  width = 5,
  height = 5,
  dpi = 1200
)  
