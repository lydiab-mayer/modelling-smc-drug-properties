library(RColorBrewer)
library(tidyr)
library(hetGP)

train_GP_E5_noninf <- function(setting) {
  
  
  EIRs = c(5,9,20,47,150)
  
  
  
  dataplotlist <- list()
  
  
  for (i in 1:length(EIRs)) {
    
    
    EIR <- EIRs[i]
    
    #import dataset
    
    save_id <- paste0("non_inferiority_",setting[, "seasonality"],"_",   EIR  ,"_",setting[, "SMC_HL"],"_",setting[, "LAI_dec"])
    seeds <- read.table(paste('/scicore/home/penny/GROUP/smc_lai/E5_Comparison/postprocessing/seeds_E5_comparison_',setting[, "seasonality"],"_",
                              EIR  ,'_',setting[, "SMC_HL"],"_",setting[, "LAI_dec"],'.txt', sep = ""), 
                        header = T, as.is = TRUE, stringsAsFactors = FALSE)
    mean(seeds$prev_beg)
    
    load("/scicore/home/penny/burlyd00/smc_lai/param_ranges_trial.RData")
    prev <- mean(seeds$prev_beg)
    
    non_inferiority_threshold <- mean(seeds$UL)
    
    #split into test and train
    n_seeds <- length(unique(seeds$seed))
    n_points = round(nrow(seeds)/n_seeds*0.80)
    
    seeds$id <- rep(seq(1,nrow(seeds)/n_seeds),each=n_seeds)
    index_train = sample(nrow(seeds)/n_seeds, n_points)
    index_test = setdiff(1:nrow(seeds)/n_seeds, index_train)
    
    train_data <-seeds[seeds$id %in% index_train,]
    test_data <- seeds[seeds$id %in% index_test,]
    test_data2 =   test_data %>% group_by(Scenario_Name) %>% summarise_at(c(names(test_data)[which(names(test_data)=="Halflife"):length(names(test_data) ) ]),mean,na.rm=TRUE)
    
    # train GP on upper limit of the confidence interval
    predicted= "CI_high_HR"
    
    train_model_matern <- train_GP_matern(train_data[ ,c(rownames(param_ranges), predicted)])
    test_model_matern  <- test_GP(GP_model=train_model_matern, train_data=train_data[,c(rownames(param_ranges), predicted)],
                                  test_data=test_data2[,c(rownames(param_ranges), predicted)])
    
    
    
    test_GP_plot(GP_model=train_model_matern, test=test_data2[,c(rownames(param_ranges), predicted)],save=paste0("/scicore/home/penny/GROUP/smc_lai/E5_Comparison/gp/",save_id,"_",predicted,".jpeg"))
    
    R2 <- cor(test_model_matern$test_data$CI_high_HR,test_model_matern$test_data$predicted_prev_red)
    
    #store output
    cv_result <- list(GP_model = train_model_matern)
    cv_result$prev <- prev
    cv_result$predicted <- predicted
    cv_result$R2 <- R2
    cv_result$train_data <- train_data
    cv_result$test_data <- test_data
    cv_result$UL =  mean(seeds$UL)
    
    cv_result$seasonality =  setting[i, "seasonality"]
    cv_result$SMC_HL =  setting[i, "SMC_HL"]
    cv_result$LAI_dec =  setting[i, "LAI_dec"]
    
    save(cv_result, file = paste0("/scicore/home/penny/GROUP/smc_lai/E5_Comparison/gp/",save_id,"_",predicted,".RData"))
    
  }}

train_GP_E5_incred <- function(setting) {
  
  
    EIRs = c(5,9,20,47,150)
  
  
  
  dataplotlist <- list()
  
  
  for (i in 1:length(EIRs)) {
    
    
    EIR <- EIRs[i]
    
    #import dataset
    
    save_id <- paste0("non_inferiority_",setting[, "seasonality"],"_",   EIR  ,"_",setting[, "SMC_HL"],"_",setting[, "LAI_dec"])
    seeds <- read.table(paste('/scicore/home/penny/GROUP/smc_lai/E5_Comparison/postprocessing/seeds_E5_comparison_',setting[, "seasonality"],"_",
                              EIR  ,'_',setting[, "SMC_HL"],"_",setting[, "LAI_dec"],'.txt', sep = ""), 
                        header = T, as.is = TRUE, stringsAsFactors = FALSE)
    
    prev <- mean(seeds$prev_beg)
    load("/scicore/home/penny/burlyd00/smc_lai/param_ranges_trial.RData")
    
    
    
    #split into test and train
    n_seeds <- length(unique(seeds$seed))
    n_points = round(nrow(seeds)/n_seeds*0.80)
    
    seeds$id <- rep(seq(1,nrow(seeds)/n_seeds),each=n_seeds)
    index_train = sample(nrow(seeds)/n_seeds, n_points)
    index_test = setdiff(1:nrow(seeds)/n_seeds, index_train)
    
    train_data <-seeds[seeds$id %in% index_train,]
    test_data <- seeds[seeds$id %in% index_test,]
    test_data2 =   test_data %>% group_by(Scenario_Name) %>% summarise_at(c(names(test_data)[which(names(test_data)=="Halflife"):length(names(test_data) ) ]),mean,na.rm=TRUE)
    
    # train GP on upper limit of the confidence interval
    predicted= "incred"
    
    train_model_matern <- train_GP_matern(train_data[ ,c(rownames(param_ranges), predicted)])
    test_model_matern  <- test_GP(GP_model=train_model_matern, train_data=train_data[,c(rownames(param_ranges), predicted)],
                                  test_data=test_data2[,c(rownames(param_ranges), predicted)])
    
    
    
    test_GP_plot(GP_model=train_model_matern, test=test_data2[,c(rownames(param_ranges), predicted)],save=paste0("/scicore/home/penny/GROUP/smc_lai/E5_Comparison/gp/",save_id,"_",predicted,".jpeg"))
    
    R2 <- cor(test_model_matern$test_data$incred,test_model_matern$test_data$predicted_prev_red)
    
    #store output
    cv_result <- list(GP_model = train_model_matern)
    cv_result$incred_SMC <- mean(seeds$incred_SMC)
    cv_result$prev <- prev
    
    cv_result$predicted <- predicted
    cv_result$R2 <- R2
    cv_result$train_data <- train_data
    cv_result$test_data <- test_data
    cv_result$UL =  mean(seeds$UL)
    
    cv_result$seasonality =  setting[i, "seasonality"]
    cv_result$SMC_HL =  setting[i, "SMC_HL"]
    cv_result$LAI_dec =  setting[i, "LAI_dec"]
    
    save(cv_result, file = paste0("/scicore/home/penny/GROUP/smc_lai/E5_Comparison/gp/",save_id,"_",predicted,".RData"))
    
  }
  
  }

extract_data <- function(cont_data) {
  data_list <- list()
  for( i in 1: length(cont_data$data)) {
    
    if(length(cont_data$data[[i]])>0) {
      df <- cont_data$data[[i]] %>%
        group_by(x) %>%
        dplyr::summarize(Halflife = mean(y, na.rm=TRUE))
      
      efficacy_val <- seq(0.7,1,by=0.01)
      
      df <- df[which(df$x %in%efficacy_val ), ]
      
      missing <- efficacy_val[-which(efficacy_val %in%df$x ) ]
      if(length(missing) >0){ 
        add_points <- data.frame(cbind(missing,160))
        names(add_points) <- c("x","Halflife")
        df <- data.frame(rbind(df,add_points))
      }
      
      df$EIR_group <- i
    } else { df=NULL }
    data_list[[i]] <- df
  }
  return(data_list)
} 

transform_data <- function(predplotlist){
  
  cont_data <-   ggplot_build(ggplot()+
                                stat_contour(data = predplotlist[[1]],aes(x= efficacy, y= halflife,z=non_inferiority, color="EIR 5"),size=2,breaks=) +
                                stat_contour(data = predplotlist[[2]],aes(x= efficacy, y= halflife,z=non_inferiority, color="EIR 20"),size=2,breaks=1 ) +
                                stat_contour(data = predplotlist[[3]],aes(x= efficacy, y= halflife,z=non_inferiority, color="EIR 70"),size=2,breaks=1 ) +
                                stat_contour(data = predplotlist[[4]],aes(x= efficacy, y= halflife,z=non_inferiority, color="EIR 100"),size=2,breaks=1 ) +
                                stat_contour(data = predplotlist[[5]],aes(x= efficacy, y= halflife,z=non_inferiority, color="EIR 200"),size=2,breaks=1 ) )
  
  
  
  data <- extract_data(cont_data)
  
  data <- bind_rows(data)
  if(length(data>0)) {
    data <- data.frame(spread(data, EIR_group, Halflife))
    
    data$X6 <- 160}
  
  return(data) 
}

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
 
   cont_data <-   ggplot_build(ggplot()+
                                stat_contour(data = predplotlist[[1]],aes(x= efficacy, y= halflife,z=incred, color="EIR 5"),size=2,breaks=breaks) +
                                stat_contour(data = predplotlist[[2]],aes(x= efficacy, y= halflife,z=incred, color="EIR 20"),size=2,breaks=breaks ) +
                                stat_contour(data = predplotlist[[3]],aes(x= efficacy, y= halflife,z=incred, color="EIR 70"),size=2,breaks=breaks ) +
                                stat_contour(data = predplotlist[[4]],aes(x= efficacy, y= halflife,z=incred, color="EIR 100"),size=2,breaks=breaks ) +
                                stat_contour(data = predplotlist[[5]],aes(x= efficacy, y= halflife,z=incred, color="EIR 200"),size=2,breaks=breaks ))
  
  
  
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
  load("~/smc_lai/param_ranges_trial.RData")
  
  
  EIRs = c(5,9,20,47,150)
  
  
  
  
  predplotlist <- list()
  dataplotlist <- list()
  for (i in 1:length(EIRs)) {
    
    
    EIR <- EIRs[i]
    
    
    
    gp_id  <- paste0("non_inferiority_",setting[, "seasonality"],"_",   EIR  ,"_",setting[, "SMC_HL"],"_",setting[, "LAI_dec"])
    
    load( file = paste0("/scicore/home/penny/GROUP/smc_lai/E5_Comparison/gp/",gp_id,"_","CI_high_HR",".RData"))
    
    
    Efficacy <- seq(param_ranges["Efficacy",1],param_ranges["Efficacy",2],0.01)
    Halflife <- seq(param_ranges["Halflife",1],param_ranges["Halflife",2],0.1)
    Xcand <- expand.grid(Halflife,Efficacy)
    names(Xcand) <- c("Halflife","Efficacy")
    
    
    prediction_for_plot <- predict(x = as.matrix(Xcand), 
                                   object = cv_result$GP_model)
    upper_limit_CI   <- data.frame(cbind(Xcand,pred=prediction_for_plot$mean ))
    
    
    
    colnames(upper_limit_CI) <- c("halflife", "efficacy", "upper_bound")
    df <- upper_limit_CI
    df$UL <- cv_result$UL
    df$non_inferiority <-   ifelse(df$upper_bound>df$UL,0, 1)
    
    predplotlist[[i]] <- df
    
    
    # load dataplots to compare predictions
    
    agg <- read.table(paste('/scicore/home/penny/GROUP/smc_lai/E5_Comparison/postprocessing/seeds_E5_comparison_',setting[, "seasonality"],"_",
                            EIR  ,'_',setting[, "SMC_HL"],"_",setting[, "LAI_dec"],'.txt', sep = ""), 
                      header = T, as.is = TRUE, stringsAsFactors = FALSE)
    q1 <-   ggplot()+
      geom_point(data= agg, aes(x= Efficacy, y= Halflife,colour=as.factor(non_inferiority)) )+
      labs( x ="Efficacy [%]", y = "Halflife [days]", colour = "")+
      #  scale_color_manual(values=    c("#999999", "#56B4E9"), labels = c("inferior", "non-inferior")) +
      theme(title = element_text(size=20),
            axis.title = element_text(size = 20),
            axis.text = element_text(size = 20),
            legend.text = element_text(size = 20))+
      
      theme(panel.border = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.text=element_text(size=20),
            axis.line = element_line(colour = "black"),
            axis.text=element_text(size=20,colour="black"),
            axis.text.x=element_text(margin = margin(t = 5)),
            axis.title=element_text(size=20,face="bold"),  
            legend.title=element_text(size=20))+
      theme(legend.background=element_blank()) + theme(legend.key=element_blank(),
                                                       legend.text = element_text(size=20,face="bold"))
    
    
    dataplotlist[[i]] <- q1
    
  }
  
  
  
  q <- ggarrange(plotlist = dataplotlist,nrow=1,common.legend=TRUE,legend="bottom")
  

  brewer.pal( n=11, name = "Blues")
  colors <- c(  "5"="#f3cad2", "20"="#e89bb9" ,"70"="#cb6ca2","100" ="#9e4387","200"= "#6d1c68")

 if(setting[, "seasonality"]=="Mali"){ 
   labels <-     c(0.4,1,1.6,2.2,2.9)
 }else{
   labels <-  c(0.5, 1, 1.5, 2.0, 2.8) 
 }
      

  
  
  data_wide <- transform_data(predplotlist)
  data_wide$x <- data_wide$x *100
  regimen = ifelse(setting[, "seasonality"]=="Mali","long season","short season")
  decay = ifelse(setting[, "LAI_dec"]=="exp","exponential LAIs",
                 ifelse( setting[, "LAI_dec"]=="hill","sigmoidal LAIs","bi-phasic LAIs"))
  
  
  p <-   ggplot()+
    {if(length(data_wide$X2)>0 ) geom_ribbon(data= data_wide, aes(x=x , ymin = X1, ymax =X2, fill = "5"))}  +
    {if(length(data_wide$X3)>0 ) geom_ribbon(data= data_wide, aes(x=x , ymin = X2, ymax =X3, fill = "20"))} +
   {if(length(data_wide$X4)>0 ) geom_ribbon(data= data_wide, aes(x=x , ymin = X3, ymax =X4, fill = "70"))} +
    {if(length(data_wide$X5)>0 )  geom_ribbon(data= data_wide, aes(x=x , ymin = X4, ymax =X5, fill = "100"))} +
    {if(length(data_wide$X5)>0 )  geom_ribbon(data= data_wide, aes(x=x , ymin = X5, ymax =X6, fill = "200"))} +
    ggtitle(paste0(regimen,"; ",decay)) +
    coord_cartesian(xlim=c(70,100),
                    ylim=c(30,150),expand=FALSE)+
    scale_fill_manual(values=colors,
                      labels=labels,
                      breaks=c("5","20","70","100","200"),drop = FALSE)+
    labs( fill=  expression(paste("initial cases per \nperson per year"["0.25-5y"])),  x ="Initial protective efficacy [%]", y = "Half-life of protective efficacy [days]")+
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
    theme(legend.background=element_blank()) + theme(legend.position = "bottom" ,legend.text = element_text(size=20,face="bold"))+guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    theme(legend.key.width = unit(2.5,"cm"))+
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
  
  
  
  
 # "1_0.8"="#faddc3", "2_0.4"="#f9c29c" ,"3_0.8"="#f6a173","4_0.8" ="#f17b51","5_0.8"= "#ea4c3b")

  
  colors <- c(  "1_0.4"="#d1fdb4", "2_0.4"="#8fccb4" ,"3_0.4"="#579c97","4_0.4" ="#2a6d7a","5_0.4"= "#0e3f5c", 
                "1_0.6"="#f3cad2", "2_0.6"="#e89bb9" ,"3_0.6"="#cb6ca2","4_0.6" ="#9e4387","5_0.6"= "#6d1c68",
                "1_0.8"="#fef5eb", "2_0.8"="#ffcda1" ,"3_0.8"="#f7921e","4_0.8" ="#d35001","5_0.8"= "#802a07")
  

  colors1 <- c(  "1"="#d1fdb4", "2"="#8fccb4" ,"3"="#579c97","4" ="#2a6d7a","5"= "#0e3f5c")
  colors2 <- c(  "1"="#f3cad2", "2"="#e89bb9" ,"3"="#cb6ca2","4" ="#9e4387","5"= "#6d1c68")
  colors3 <- c( "1"="#fef5eb", "2"="#ffcda1" ,"3"="#f7921e","4" ="#d35001","5"= "#802a07")
  
  
    data_wide <- transform_data_incred(predplotlist)

  data_wide$cols <- paste0(data_wide$EIR_group,"_",data_wide$level)
   if(setting[, "seasonality"]=="Mali"){ 
    labels <-     c(0.4,1,1.6,2.2,2.9)
  }else{
    labels <-  c(0.5, 1, 1.5, 2.0, 2.8) 
  }
  
  regimen = ifelse(setting[, "seasonality"]=="Mali","long season","short season")
  decay = ifelse(setting[, "LAI_dec"]=="exp","exponential LAIs",
                 ifelse( setting[, "LAI_dec"]=="hill","sigmoidal LAIs","bi-phasic LAIs"))
  
  
  
  ###test 
  
  
 
  df_leg <- data.frame(cbind("40 %"=labels,"60 %"=labels,"80 %"=labels ))
rownames(df_leg) <- labels
colnames(df_leg) <- c("40 %", " 60 %", "80 %")



a <- melt(df_leg)
a$fill <- names(colors)

a$valuefac <- as.factor(a$value)
 legend <-  ggplot(data= a, aes(variable,valuefac)) +
  geom_tile(aes(fill = fill)) +
  # take out dt$...
  # take out 3 from aes
  scale_fill_manual(values = c(colors)) +
   labs(  
          x ="Estimated incidence reduction [%]", y =  expression(bold(paste("Initial cases per person per year"["0.25-5y"]))))+
    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.text=element_text(size=15),
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=25,colour="black", face = "bold"),
          axis.text.x=element_text(margin = margin(t = 5)),
          axis.title=element_text(size=25,face="bold"), 
          legend.title=element_text(size=25,face="bold"))+ theme(legend.key=element_blank()) +
    theme(legend.background=element_blank()) + theme(legend.text = element_text(size=15,face="bold"))+
    theme(legend.key.width = unit(4,"cm")) + theme(legend.key=element_blank(),  legend.text = element_text(size=20,face="bold"))+
    theme(plot.title = element_text(size = 25, face = "bold",hjust = 0.2, vjust=2.12))+   
    theme(legend.position =  "none")

  
  
 
  
  #### test ende
  leg1 <- get_legend(ggplot()+ geom_line(data= data_wide[which(data_wide$level==0.6),], aes(x=x*100 , y = Halflife, 
                                                                                            color = as.character(EIR_group)),show.legend = TRUE,
                                         size=2 ) +
                       labs( color =   expression(paste("40% incidence \nreduction","\n initial cpppy"["0.25-5y"], sep="")), 
                             x ="Initial protective efficacy [%]", y = "Half-life of protective efficacy [days]")+
                       scale_color_manual(values=colors1,
                                          labels=labels,
                                          breaks=c("1","2","3","4","5"),drop = FALSE)+
                       theme( legend.text=element_text(size=20),
                              legend.title=element_text(size=25,face="bold"))+ theme(legend.key=element_blank(), legend.position = "right") +
                       theme(legend.background=element_blank()) + theme(legend.text = element_text(size=20,face="bold"))+
                       theme(legend.key.width = unit(4,"cm"))+ theme(legend.key=element_blank(),  legend.text = element_text(size=20,face="bold")))
                      
  
 
  
  leg2 <- get_legend(ggplot()+ geom_line(data= data_wide[which(data_wide$level==0.6),], aes(x=x*100 , y = Halflife, 
                                                                                            color = as.character(EIR_group)),show.legend = TRUE,
                                         size=2 ) +
                       labs( color =   expression(paste("60% incidence \nreduction","initial cpppy"["0.25-5y"], sep="")), 
                             x ="Initial protective efficacy [%]", y = "Half-life of protective efficacy [days]")+
                       scale_color_manual(values=colors2,
                                          labels=labels,
                                          breaks=c("1","2","3","4","5"),drop = FALSE)+
                       theme( legend.text=element_text(size=20),
                              legend.title=element_text(size=25,face="bold"))+ theme(legend.key=element_blank(), legend.position = "right") +
                       theme(legend.background=element_blank()) + theme(legend.text = element_text(size=20,face="bold"))+
                       theme(legend.key.width = unit(4,"cm"))+ theme(legend.key=element_blank(),  legend.text = element_text(size=20,face="bold")))
  leg3 <- get_legend(ggplot()+ geom_line(data= data_wide[which(data_wide$level==0.6),], aes(x=x*100 , y = Halflife, 
                                                                                            color = as.character(EIR_group)),show.legend = TRUE,
                                         size=2 ) +
                       labs( color =   expression(paste("80% incidence \nreduction","initial cpppy"["0.25-5y"], sep="")), 
                             x ="Initial protective efficacy [%]", y = "Half-life of protective efficacy [days]")+
                       scale_color_manual(values=colors3,
                                          labels=labels,
                                          breaks=c("1","2","3","4","5"),drop = FALSE)+
                       theme( legend.text=element_text(size=20),
                              legend.title=element_text(size=25,face="bold"))+ theme(legend.key=element_blank(), legend.position = "right") +
                       theme(legend.background=element_blank()) + theme(legend.text = element_text(size=20,face="bold"))+
                       theme(legend.key.width = unit(4,"cm"))+ theme(legend.key=element_blank(),  legend.text = element_text(size=20,face="bold")))
  
    
  
  
  p <-   
   # {if(length(data_wide$X1)>0 ) geom_line(data= data_wide, aes(x=x , y = X1,  color = "5", linetype=incred),size=1.5,)}  +
    #{if(length(data_wide$X2)>0 ) geom_line(data= data_wide, aes(x=x , y = X2,  color = "20", linetype=incred),size=1.5,)}  +
   #{if(length(data_wide$X3)>0 ) geom_line(data= data_wide, aes(x=x , y = X3,  color = "70", linetype=incred),size=1.5,)}  +
    #{if(length(data_wide$X4)>0 ) geom_line(data= data_wide, aes(x=x , y = X4,  color = "100", linetype=incred),size=1.5,)}  +
    #{if(length(data_wide$X5)>0 ) geom_line(data= data_wide, aes(x=x , y = X5,  color = "200", linetype=incred),size=1.5)}  +
    ggplot()+ geom_line(data= data_wide, aes(x=x*100 , y = Halflife, 
                                             color = as.character(cols)) ,size=2)  +
  
    coord_cartesian(xlim=c(70,100),
                    ylim=c(30,150),expand=FALSE)+
    ggtitle(paste0(regimen,"; ",decay)) +
    scale_color_manual(values=colors,
                      labels=rep(labels,3),
                      breaks=rep(c("1","2","3","4","5"),3),drop = FALSE)+
    labs( color =   expression(paste("initial cpppy"["0.25-5y"])), 
           x ="Initial protective efficacy [%]", y = "Half-life of protective efficacy [days]")+
     theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.text=element_text(size=15),
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=25,colour="black"),
          axis.text.x=element_text(margin = margin(t = 5)),
          axis.title=element_text(size=25,face="bold"), 
          legend.title=element_text(size=20,face="bold"))+ theme(legend.key=element_blank()) +
    theme(legend.background=element_blank()) + theme(legend.text = element_text(size=15,face="bold"))+
    theme(legend.key.width = unit(4,"cm")) + theme(legend.key=element_blank(),  legend.text = element_text(size=20,face="bold"))+
    theme(plot.title = element_text(size = 25, face = "bold",hjust = 0.2, vjust=2.12))
  
  

  # combine legend 1 & 2
  leg12 <- plot_grid(leg1, leg2,leg3,ncol=3
  )
 
  
  ggsave(
   filename= paste0("./Outputs/E5_Compare/incredcontours_",gp_id,".jpg"),
    plot = last_plot(),

    width = 7,
    height = 5,
    dpi = 600,
  
  )
  return(list(p,legend) )
  
  
}


### plots for thesis presentation 

pred_GP_E5_noninf_talk <- function(setting) {
  load("~/smc_lai/param_ranges_trial.RData")
  
  
  EIRs = c(5,9,20,47,150)
  
  
  
  
  predplotlist <- list()
  dataplotlist <- list()
  for (i in 1:length(EIRs)) {
    
    
    EIR <- EIRs[i]
    
    
    
    gp_id  <- paste0("non_inferiority_",setting[, "seasonality"],"_",   EIR  ,"_",setting[, "SMC_HL"],"_",setting[, "LAI_dec"])
    
    load( file = paste0("/scicore/home/penny/GROUP/smc_lai/E5_Comparison/gp/",gp_id,"_","CI_high_HR",".RData"))
    
    
    Efficacy <- seq(param_ranges["Efficacy",1],param_ranges["Efficacy",2],0.01)
    Halflife <- seq(param_ranges["Halflife",1],param_ranges["Halflife",2],0.1)
    Xcand <- expand.grid(Halflife,Efficacy)
    names(Xcand) <- c("Halflife","Efficacy")
    
    
    prediction_for_plot <- predict(x = as.matrix(Xcand), 
                                   object = cv_result$GP_model)
    upper_limit_CI   <- data.frame(cbind(Xcand,pred=prediction_for_plot$mean ))
    
    
    
    colnames(upper_limit_CI) <- c("halflife", "efficacy", "upper_bound")
    df <- upper_limit_CI
    df$UL <- cv_result$UL
    df$non_inferiority <-   ifelse(df$upper_bound>df$UL,0, 1)
    
    predplotlist[[i]] <- df
    
    
    # load dataplots to compare predictions
    
    agg <- read.table(paste('/scicore/home/penny/GROUP/smc_lai/E5_Comparison/postprocessing/seeds_E5_comparison_',setting[, "seasonality"],"_",
                            EIR  ,'_',setting[, "SMC_HL"],"_",setting[, "LAI_dec"],'.txt', sep = ""), 
                      header = T, as.is = TRUE, stringsAsFactors = FALSE)
    q1 <-   ggplot()+
      geom_point(data= agg, aes(x= Efficacy, y= Halflife,colour=as.factor(non_inferiority)) )+
      labs( x ="Efficacy [%]", y = "Halflife [days]", colour = "")+
      #  scale_color_manual(values=    c("#999999", "#56B4E9"), labels = c("inferior", "non-inferior")) +
      theme(title = element_text(size=20),
            axis.title = element_text(size = 20),
            axis.text = element_text(size = 20),
            legend.text = element_text(size = 20))+
      
      theme(panel.border = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.text=element_text(size=20),
            axis.line = element_line(colour = "black"),
            axis.text=element_text(size=20,colour="black"),
            axis.text.x=element_text(margin = margin(t = 5)),
            axis.title=element_text(size=20,face="bold"),  
            legend.title=element_text(size=20))+
      theme(legend.background=element_blank()) + theme(legend.key=element_blank(),
                                                       legend.text = element_text(size=20,face="bold"))
    
    
    dataplotlist[[i]] <- q1
    
  }
  
  
  
  q <- ggarrange(plotlist = dataplotlist,nrow=1,common.legend=TRUE,legend="bottom")
  
  
  brewer.pal( n=9, name = "Blues")
  colors <- c(  "70"="#cb6ca2")
  
    labels <-     c(1)
 
  
  
  data_wide <- transform_data(predplotlist)
  data_wide$x <- data_wide$x *100
  regimen = ifelse(setting[, "seasonality"]=="Mali","long season","short season")
  decay = ifelse(setting[, "LAI_dec"]=="exp","exponential LAIs",
                 ifelse( setting[, "LAI_dec"]=="hill","sigmoidal LAIs","bi-phasic LAIs"))
  
  
  p <-   ggplot()+
  #  {if(length(data_wide$X2)>0 ) geom_ribbon(data= data_wide, aes(x=x , ymin = X1, ymax =X2, fill = "5"))}  +
    {if(length(data_wide$X3)>0 ) geom_ribbon(data= data_wide, aes(x=x , ymin = X2, ymax =X6, fill = "70"))} +
    #{if(length(data_wide$X4)>0 ) geom_ribbon(data= data_wide, aes(x=x , ymin = X3, ymax =X6, fill = "70"))} +
    #{if(length(data_wide$X5)>0 )  geom_ribbon(data= data_wide, aes(x=x , ymin = X4, ymax =X5, fill = "100"))} +
    #{if(length(data_wide$X5)>0 )  geom_ribbon(data= data_wide, aes(x=x , ymin = X5, ymax =X6, fill = "200"))} +
    ggtitle(paste0(regimen,"; ",decay)) +
    coord_cartesian(xlim=c(70,100),
                    ylim=c(30,150),expand=FALSE)+
    scale_fill_manual(values=colors,
                      labels=labels,
                      breaks=c("70"),drop = FALSE)+
    labs( fill=  expression(paste("initial cpppy"["0.25-5y"])),  x ="Initial protective efficacy [%]", y = "Half-life of protective efficacy [days]")+
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
    theme(legend.background=element_blank()) + theme(legend.position = "bottom" ,legend.text = element_text(size=20,face="bold"))+guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    theme(legend.key.width = unit(2.5,"cm"))+
    theme(plot.title = element_text(size = 25, face = "bold",hjust = 0.2, vjust=2.12))
  
  
  

  
  return(p)
  
  
}

pred_GP_E5_incred_talk <- function(setting) {
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
  
  
  
  
  # "1_0.8"="#faddc3", "2_0.4"="#f9c29c" ,"3_0.8"="#f6a173","4_0.8" ="#f17b51","5_0.8"= "#ea4c3b")
  
  
  
  colors <- c(  "3_0.4"="#579c97",
               "3_0.6"="#cb6ca2",
               "3_0.8"="#f7921e")
  
  
  colors1 <- c(  "3"="#579c97")
  colors2 <- c( "3"="#cb6ca2")
  colors3 <- c("3"="#f7921e")
  

  
  
  
  data_wide <- transform_data_incred(predplotlist)
  
  data_wide$cols <- paste0(data_wide$EIR_group,"_",data_wide$level)
 
  data_wide <- subset(data_wide, EIR_group==3)
    labels <-     c(1.6)
  
  regimen = ifelse(setting[, "seasonality"]=="Mali","long season","short season")
  decay = ifelse(setting[, "LAI_dec"]=="exp","exponential LAIs",
                 ifelse( setting[, "LAI_dec"]=="hill","sigmoidal LAIs","bi-phasic LAIs"))
  leg1 <- get_legend(ggplot()+ geom_line(data= data_wide[which(data_wide$level==0.6),], aes(x=x*100 , y = Halflife, 
                                                                                            color = as.character(EIR_group)),show.legend = TRUE,
                                         size=2 ) +
                       labs( color =   expression(paste("40% incidence\nreduction:","initial cpppy"["0.25-5y "] , sep="")), 
                             x ="Initial protective efficacy [%]", y = "Half-life of protective efficacy [days]")+
                       scale_color_manual(values=colors1,
                                          labels=labels,
                                          breaks=c("3"),drop = FALSE)+
                       theme( legend.text=element_text(size=15),
                              legend.title=element_text(size=20,face="bold"))+ theme(legend.key=element_blank(), legend.position = "right") +
                       theme(legend.background=element_blank()) + theme(legend.text = element_text(size=15,face="bold"))+
                       theme(legend.key.width = unit(4,"cm"))+ theme(legend.key=element_blank(),  legend.text = element_text(size=20,face="bold")))
  
  
  leg2 <- get_legend(ggplot()+ geom_line(data= data_wide[which(data_wide$level==0.6),], aes(x=x*100 , y = Halflife, 
                                                                                            color = as.character(EIR_group)),show.legend = TRUE,
                                         size=2 ) +
                       labs( color =   expression(paste("60% incidence\nreduction:","initial cpppy"["0.25-5y "] , sep="")), 
                             x ="Initial protective efficacy [%]", y = "Half-life of protective efficacy [days]")+
                       scale_color_manual(values=colors2,
                                          labels=labels,
                                          breaks=c("3"),drop = FALSE)+
                       theme( legend.text=element_text(size=15),
                              legend.title=element_text(size=20,face="bold"))+ theme(legend.key=element_blank(), legend.position = "right") +
                       theme(legend.background=element_blank()) + theme(legend.text = element_text(size=15,face="bold"))+
                       theme(legend.key.width = unit(4,"cm"))+ theme(legend.key=element_blank(),  legend.text = element_text(size=20,face="bold")))
  leg3 <- get_legend(ggplot()+ geom_line(data= data_wide[which(data_wide$level==0.6),], aes(x=x*100 , y = Halflife, 
                                                                                            color = as.character(EIR_group)),show.legend = TRUE,
                                         size=2 ) +
                       labs( color =   expression(paste("80% incidence\nreduction:","initial cpppy"["0.25-5y "] , sep="")), 
                             x ="Initial protective efficacy [%]", y = "Half-life of protective efficacy [days]")+
                       scale_color_manual(values=colors3,
                                          labels=labels,
                                          breaks=c("3"),drop = FALSE)+
                       theme( legend.text=element_text(size=15),
                              legend.title=element_text(size=20,face="bold"))+ theme(legend.key=element_blank(), legend.position = "right") +
                       theme(legend.background=element_blank()) + theme(legend.text = element_text(size=15,face="bold"))+
                       theme(legend.key.width = unit(4,"cm"))+ theme(legend.key=element_blank(),  legend.text = element_text(size=20,face="bold")))
  
  
  
  
  p <-   
    # {if(length(data_wide$X1)>0 ) geom_line(data= data_wide, aes(x=x , y = X1,  color = "5", linetype=incred),size=1.5,)}  +
    #{if(length(data_wide$X2)>0 ) geom_line(data= data_wide, aes(x=x , y = X2,  color = "20", linetype=incred),size=1.5,)}  +
    #{if(length(data_wide$X3)>0 ) geom_line(data= data_wide, aes(x=x , y = X3,  color = "70", linetype=incred),size=1.5,)}  +
    #{if(length(data_wide$X4)>0 ) geom_line(data= data_wide, aes(x=x , y = X4,  color = "100", linetype=incred),size=1.5,)}  +
    #{if(length(data_wide$X5)>0 ) geom_line(data= data_wide, aes(x=x , y = X5,  color = "200", linetype=incred),size=1.5)}  +
    ggplot()+ geom_line(data= data_wide, aes(x=x*100 , y = Halflife, 
                                             color = as.character(cols)) ,size=2)  +
    
    coord_cartesian(xlim=c(70,100),
                    ylim=c(30,150),expand=FALSE)+
    ggtitle(paste0(regimen,"; ",decay)) +
    scale_color_manual(values=colors,
                       labels=rep(labels,3),
                       breaks=c("2","2","2"),drop = FALSE)+
    labs( color =   expression(paste("initial cpppy"["0.25-5y"])), 
          x ="Initial protective efficacy [%]", y = "Half-life of protective efficacy [days]")+
    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.text=element_text(size=15),
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=20,colour="black"),
          axis.text.x=element_text(margin = margin(t = 5)),
          axis.title=element_text(size=20,face="bold"), 
          legend.title=element_text(size=20,face="bold"))+ theme(legend.key=element_blank()) +
    theme(legend.background=element_blank()) + theme(legend.text = element_text(size=15,face="bold"))+
    theme(legend.key.width = unit(4,"cm"))+ theme(legend.key=element_blank(),  legend.text = element_text(size=20,face="bold"))+
    theme(plot.title = element_text(size = 25, face = "bold",hjust = 0.2, vjust=2.12))
  
  
  
  # combine legend 1 & 2
  leg12 <- plot_grid(leg1, leg2,leg3,nrow=3
  )
  
  
 
  return(list(p,leg12) )
  
  
}
