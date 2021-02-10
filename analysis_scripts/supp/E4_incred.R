library(RColorBrewer)
library(tidyr)
library(hetGP)



extract_data_incred <- function(cont_data) {
  data_list <- list()
  for( i in 1: length(cont_data$data)) {
    
    if(length(cont_data$data[[i]])>0) {
      efficacy_val <- seq(0.7,1,by=0.01)
      df <-subset(cont_data$data[[i]], x%in%efficacy_val )
      
      df <- df %>%
        group_by(x,level) %>%
        dplyr::summarize(Halflife = mean(y, na.rm=TRUE))
      

add_missingpoints <- function(x) {
  efficacy_val <- seq(0.7,1,by=0.01)

      x <- x[which(x$x %in%efficacy_val ),c("x","Halflife") ]
      
      missing <- efficacy_val[-which(efficacy_val %in%x$x ) ]
      if(length(missing) >0){ 
        add_points <- data.frame(cbind(missing,160))
        names(add_points) <- c("x","Halflife")
        x <- data.frame(rbind(x,add_points))

      } else {x <- x}
      return(x)
    }
    
df_list <- df %>%
  group_by(level) %>%
  group_walk(~add_missingpoints(.x))

df$EIR_group <- i
    } else { df=NULL }
    data_list[[i]] <- df
  }
  return(data_list)
} 


transform_data_incred<- function(predplotlist){
  breaks <- c(0.4,0.6,0.8)
 
  for (i in length(predplotlist)) {
    df <- predplotlist[[i]]
    
 dfs <-    split(df, with(df, Coverage), drop = TRUE)
 
 cont_data <-   ggplot_build(
   ggplot()+
     stat_contour(data = dfs[[1]],aes(x= Efficacy, y= Halflife,z=incred, color="EIR 5"),size=2,breaks=breaks) +
     stat_contour(data = dfs[[2]],aes(x= Efficacy, y= Halflife,z=incred, color="EIR 20"),size=2,breaks=breaks ) +
     stat_contour(data = dfs[[3]],aes(x= Efficacy, y= Halflife,z=incred, color="EIR 70"),size=2,breaks=breaks ) +
     stat_contour(data = dfs[[4]],aes(x= Efficacy, y= Halflife,z=incred, color="EIR 100"),size=2,breaks=breaks ) 
 )
 
 
  }
  
   
  
  
  data <- extract_data_incred(cont_data)
  
  
  data <- bind_rows(data)
  #if(length(data>0)) {
  #  data <- data.frame(spread(data, EIR_group,Halflife))
 # data$X6 <- 160}
  
  x <- seq(0.7,1)
  level <- c(0.4,0.6,0.8)
  comb <- expand.grid(x,level)
  names(comb) <- c("x","level")
  
missing <-   data %>%
    group_by(EIR_group,level) %>%
    group_modify(~ head(.x, 1L))
missing <- subset(missing, x>0.7)
missing$x <- missing$x-0.01
missing$Halflife <- 160
  

data <- data.frame(rbind(data,missing))
  return(data) 
}

pred_GP_E5_noninf <- function(setting) {
  load("~/smc_lai/param_ranges.RData")
  
  
points_low <-c(3,4,8,28,150)
points_high <- c(5,9,20,47,150)

 if(setting$Access==0.1) {
   EIRs <- points_low
 labels <-  c(0.45,  0.72,1.3,  2.3, 3.2) 
 title <- "low access to healthcare"
 }else{
   labels <-     c(0.4,1,1.6,2.2,2.9)
   title <- "high access to healthcare"
   EIRs <- points_high}


  
  predplotlist <- list()
  
  for (i in 1:length(EIRs)) {
    
    
    EIR <- EIRs[i]
    
    
    
    gp_id  <- paste0(setting[, "seasonality"],"_",   setting[, "IntAge"]  ,"_",setting[, "LAI_dec"],"_",setting[, "Access"],"_",EIR)
    
    load( file = paste0("/scicore/home/penny/GROUP/smc_lai/E3E4_lb/gp/pppy_LAI/gp_trained_pppy_LAI_",gp_id,".RData"))
    
    
    Efficacy <- seq(param_ranges["Efficacy",1],param_ranges["Efficacy",2],0.01)
    Halflife <- seq(param_ranges["Halflife",1],param_ranges["Halflife",2],0.1)
    Coverage <- c(0.4,0.6,0.8,1)
    Coverage_SMC <- 0.01
    
    Xcand <- expand.grid(Coverage_SMC,Efficacy,Halflife,Coverage)
    names(Xcand) <- c("Coverage_SMC","Efficacy","Halflife","Coverage")
    
    
    prediction_for_plot <- predict(x = as.matrix(Xcand), 
                                   object = GP_model)
    pred   <- data.frame(cbind(Xcand,pred=prediction_for_plot$mean ))
    
    
    
    colnames(pred) <- c("Coverage_SMC","Efficacy","Halflife","Coverage","cpppy")
   pred$incred <- (labels[i]-pred$cpppy)/labels[i]
    predplotlist[[i]] <- pred
}
  
  
  


  brewer.pal( n=9, name = "Blues")
  colors <- c(  "5"="#9ECAE1", "20"="#6BAED6" ,"70"="#4292C6","100" ="#2171B5","200"= "#08519C")
  
 if(setting[, "seasonality"]=="Mali"){ 
   labels <-     c(0.4,1,1.6,2.2,2.9)
 }else{
   labels <-  c(0.5, 1, 1.5, 2.0, 2.8) 
 }
      

  
  
  data_wide <- transform_data_incred(predplotlist)
  data_wide$x <- data_wide$x *100
  regimen = ifelse(setting[, "seasonality"]=="Mali","long season","short season")
  decay = ifelse(setting[, "LAI_dec"]=="exp","LAI 1: exponential",
                 ifelse( setting[, "LAI_dec"]=="hill","LAI 3: sigmoid","LAI 2: bi-phasic"))
  
  
  p <-   ggplot()+
    {if(length(data_wide$X2)>0 ) geom_ribbon(data= data_wide, aes(x=x , ymin = X1, ymax =X2, fill = "5"))}  +
    {if(length(data_wide$X3)>0 ) geom_ribbon(data= data_wide, aes(x=x , ymin = X2, ymax =X3, fill = "20"))} +
   {if(length(data_wide$X4)>0 ) geom_ribbon(data= data_wide, aes(x=x , ymin = X3, ymax =X4, fill = "70"))} +
    {if(length(data_wide$X5)>0 )  geom_ribbon(data= data_wide, aes(x=x , ymin = X4, ymax =X5, fill = "100"))} +
    {if(length(data_wide$X5)>0 )  geom_ribbon(data= data_wide, aes(x=x , ymin = X5, ymax =X6, fill = "200"))} +
    ggtitle(paste0(regimen,"; ", decay) )+
    coord_cartesian(xlim=c(70,100),
                    ylim=c(30,150),expand=FALSE)+
    scale_fill_manual(values=colors,
                      labels=labels,
                      breaks=c("5","20","70","100","200"),drop = FALSE)+
    labs( fill=  expression(paste("cpppy"["0.25-5y"])),  x ="Initial protective efficacy [%]", y = "Halflife of protective efficacy [d]")+
    annotate("text", x = 75, y =35 , label = "inferiority",size=10,face="bold")+
    annotate("text", x = 90, y =140 , label = "non-inferiority",size=10,face="bold")+
    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.text=element_text(size=15),
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=20,colour="black"),
          axis.text.x=element_text(margin = margin(t = 5)),
          axis.title=element_text(size=20,face="bold"), 
          legend.title=element_text(size=20,face="bold"))+ 
    theme(legend.background=element_blank()) + theme(legend.text = element_text(size=20,face="bold"))+
    theme(legend.key.width = unit(1.5,"cm"))+
   theme(plot.title = element_text(size = 25, face = "bold",hjust = 0.2, vjust=2.12))
  

  
  ggsave(
    filename =  paste0("./Outputs/E5_Compare/EIRcontours_",gp_id,".jpg"),
    plot = last_plot(),
    
    width = 7,
    height = 5,
    dpi = 600,
    
  )
  
  return(p)
  
  
}

pred_GP_E5_incred <- function(setting) {
  load("~/smc_lai/param_ranges_trial.RData")
  
  
  EIRs = c(5,9,20,47,150)
  
  
  
  
  predplotlist <- list()
  for (i in 1:length(EIRs)) {
    
    
    EIR <- EIRs[i]
    
    
    
    gp_id  <- paste0("non_inferiority_",setting[, "seasonality"],"_",   EIR  ,"_",setting[, "SMC_HL"],"_",setting[, "LAI_dec"])
    
    load( file = paste0("/scicore/home/penny/GROUP/smc_lai/E5_Comparison/gp/",gp_id,"_","incred",".RData"))
    
    
    Efficacy <- seq(param_ranges["Efficacy",1],param_ranges["Efficacy",2],0.01)
    Halflife <- seq(param_ranges["Halflife",1],param_ranges["Halflife",2],0.1)
    Xcand <- expand.grid(Halflife,Efficacy)
    names(Xcand) <- c("Halflife","Efficacy")
    
    
    prediction_for_plot <- predict(x = as.matrix(Xcand), 
                                   object = cv_result$GP_model)
    incred   <- data.frame(cbind(Xcand,pred=prediction_for_plot$mean ))
    
    
    

    colnames(incred) <- c("halflife", "efficacy", "incred")
    incred$EIR <- EIR
    incred$incred_SMC <-0.85 #cv_result$incred_SMC
    
    predplotlist[[i]] <- incred
    
    
  

  }
  
  
  
  
  brewer.pal( n=9, name = "Blues")
  colors <- c(  "1"="#C6DBEF", "2"="#6BAED6" ,"3"="#4292C6","4" ="#2171B5","5"= "#08306B")
  
  
  
  data_wide <- transform_data_incred(predplotlist)

   if(setting[, "seasonality"]=="Mali"){ 
    labels <-     c(0.4,1,1.6,2.2,2.9)
  }else{
    labels <-  c(0.5, 1, 1.5, 2.0, 2.8) 
  }
  
  regimen = ifelse(setting[, "seasonality"]=="Mali","long season","short season")
  decay = ifelse(setting[, "LAI_dec"]=="exp","LAI 1: exponential",
                 ifelse( setting[, "LAI_dec"]=="hill","LAI 3: sigmoid","LAI 2: bi-phasic"))
  p <-   
   # {if(length(data_wide$X1)>0 ) geom_line(data= data_wide, aes(x=x , y = X1,  color = "5", linetype=incred),size=1.5,)}  +
    #{if(length(data_wide$X2)>0 ) geom_line(data= data_wide, aes(x=x , y = X2,  color = "20", linetype=incred),size=1.5,)}  +
   #{if(length(data_wide$X3)>0 ) geom_line(data= data_wide, aes(x=x , y = X3,  color = "70", linetype=incred),size=1.5,)}  +
    #{if(length(data_wide$X4)>0 ) geom_line(data= data_wide, aes(x=x , y = X4,  color = "100", linetype=incred),size=1.5,)}  +
    #{if(length(data_wide$X5)>0 ) geom_line(data= data_wide, aes(x=x , y = X5,  color = "200", linetype=incred),size=1.5)}  +
    ggplot()+ geom_line(data= data_wide, aes(x=x*100 , y = Halflife, 
                                             color = as.character(EIR_group),
                                             linetype=as.factor(level)),size=1.5)  +
  
    coord_cartesian(xlim=c(70,100),
                    ylim=c(30,150),expand=FALSE)+
    scale_linetype_manual(values=c("solid","twodash", "dotted"))+
    ggtitle(paste0(regimen,"; ", decay) )+
    scale_color_manual(values=colors,
                      labels=labels,
                      breaks=c("1","2","3","4","5"),drop = FALSE)+
    labs( color =   expression(paste("cpppy"["0.25-5y"])), 
          linetype = "incidence \n reduction", x ="Initial protective efficacy [%]", y = "Halflife of protective efficacy [d]")+
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
    theme(plot.title = element_text(size = 25, face = "bold",hjust = 0.2, vjust=2.12))
  
  

  ggsave(
   filename= paste0("./Outputs/E5_Compare/incredcontours_",gp_id,".jpg"),
    plot = last_plot(),

    width = 7,
    height = 5,
    dpi = 600,
  
  )
  return(p)
  
  
}

