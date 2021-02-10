#=========================================
#
#	Analysis of OM experiments smc_lai
# -plotting cases averted in same setting example plot
# 
#	author:  Lydia Burgert
#			
# input: 
# output: 
#=========================================

library("plot3D")


res <- import_procsim(Exp="E3_SMCSMC",follow_up=5,setting="seasonal_Low_indoor",agg=TRUE,scen_id="all")

res2 <- import_procsim(Exp="E4_SMCLAI",follow_up=5,setting="seasonal_Low_indoor",agg=TRUE,scen_id="all")

# plot cases per person per year in intervention age gorup
df1 <- res[["ClinCases"]][,c("EIR","Coverage","cases","year")]
df1 <- subset(df1, year==10)
df2 <- res2[["ClinCases"]][,c("EIR","Coverage","Efficacy","Halflife","cases","year")]
df2 <- subset(df2, year==10)

df <- data.frame(cbind(df2,diff=df2$cases-df1$cases))

df$binEIR <- cut(df$EIR, breaks=c(0,5,10,15,25), labels=c("EIR 1-5","EIR 5-10","EIR 10-15","EIR 15-25"))
df$binCov <- cut(df$Coverage, breaks=c(0.0,0.25,0.50,0.75,1), labels=c("Cov 0-25","Cov 25-50","Cov 50-75","Cov 75-100"))

df <- subset(df, diff!=0)

# lai - smc (2000-1000)--> grösstenteils positiv heisst lai hat mehr fälle weil die different sonst negativ wäre--
# differenz ist nur negativ wenn: 
df$diff2 <- ifelse(df$diff>5,0,1)
ggplot(data = df, aes(x = Halflife,y=Efficacy,colour=diff2)) +
  geom_point() + facet_wrap(binEIR~binCov )+ scale_color_gradient(low = "blue", high = "red")


ggplot(data = df, aes(x = Halflife,y=Efficacy,colour=diff)) +
  geom_point() + facet_wrap(binEIR~binCov )+ scale_color_gradient(low = "blue", high = "red")

 
  
# plot cases per person per year in intervention age gorup
df1 <- res[["meanpppy"]][,c("EIR","Coverage","mean")]
df2 <- res2[["meanpppy"]][,c("EIR","Coverage","Efficacy","Halflife","mean")]

df <- data.frame(cbind(df2,diff=df2$mean-df1$mean))

df$binEIR <- cut(df$EIR, breaks=c(0,5,10,15,25), labels=c("EIR 1-5","EIR 5-10","EIR 10-15","EIR 15-25"))
df$binCov <- cut(df$Coverage, breaks=c(0.0,0.25,0.50,0.75,1), labels=c("Cov 0-25","Cov 25-50","Cov 50-75","Cov 75-100"))

df <- subset(df, diff!=0)

ggplot(data = df, aes(x = Halflife,y=Efficacy,colour=diff)) +
  geom_point() + facet_wrap(binEIR~binCov )+ scale_color_gradient(low = "blue", high = "red")

