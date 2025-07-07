
###############################################
#                  R Packages                 #
###############################################
#real dataset
library(nlme)
library(FlexParamCurve)
#plot
library(ggplot2)
library(ggbreak)
library(ggthemes)
library(fields)

#Functions
source("FunctionPopulationGrowthFDE_Conformable.R",local=FALSE)
source("FunctionPopulationGrowthODE.R",local=FALSE)
source("FunctionPopulationGrowthFDE_Caputo.R",local=FALSE)
source("FunctionNegativeLogLikFDE_Conformable.R",local=FALSE)
source("FunctionNegativeLogLikODE.R",local=FALSE)
source("FunctionNegativeLogLikFDE_Caputo.R",local=FALSE)
source("FunctionNegativeLogProfileLikFDE_Conformable.R",local=FALSE)
source("FunctionNegativeLogProfileLikFDE_Caputo.R",local=FALSE)
source("FunctionNegativeLogProfileLikFDEx0_Conformable.R",local=FALSE)
source("FunctionNegativeLogProfileLikFDEx0_Caputo.R",local=FALSE)
source("FunctionNegativeLogProfileLikFDEK_Conformable.R",local=FALSE)
source("FunctionNegativeLogProfileLikFDEK_Caputo.R",local=FALSE)
source("FunctionNegativeLogProfileLikFDElambda_Conformable.R",local=FALSE)
source("FunctionNegativeLogProfileLikFDElambda_Caputo.R",local=FALSE)
source("FunctionProfileLCI.R",local=FALSE)


###############################################
#               Caputo settings               #
###############################################
h <- 0.1
m0 <- 490
final_time <- h*m0
time <- seq(0,final_time,h)

###############################################
#               Dataset selection             #
###############################################
head(penguin.data)

#plot 1: ckage (days) versus mass (grams) : 2000 y 2002
ggplot(penguin.data,aes(x=ckage, y=weight,shape =factor(year))) +
  geom_point(size=1) +
  scale_shape_manual(values=c(1,2),name = "Year") +
  theme_classic() +
  labs(title="Field data on growth of little penguins (see Chiaradia & Nisbet 2006)",
       x="Chick age (days)", y="Body mass (grams)") +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=13,face="bold"),
        legend.title = element_text(size=13,face="bold"),
        legend.text=element_text(size=13),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major.y = element_line(),
        panel.grid.major.x = element_line())


#plot 2: Plotted values are daily means for ckage (days)
data2000<-subset(penguin.data,year=="2000",select = c("ckage","weight"))
data2002<-subset(penguin.data,year=="2002",select = c("ckage","weight"))
MeanData2000<-aggregate(data2000[,2], list(data2000$ckage), mean)
MeanData2002<-aggregate(data2002[,2], list(data2002$ckage), mean)
colnames(MeanData2000)<-c("ckage","weight")
colnames(MeanData2002)<-c("ckage","weight")
MeanData2000$year<-rep("2000",dim(MeanData2000)[1])
MeanData2002$year<-rep("2002",dim(MeanData2002)[1])
dfMean<-rbind(MeanData2000,MeanData2002)

ggplot(dfMean,aes(x=ckage, y=weight,shape =factor(year))) +
  geom_point(size=1) +
  scale_shape_manual(values=c(1,2),name = "Year") +
  theme_classic() +
  labs(title="Field data on growth of little penguins (see Chiaradia & Nisbet 2006)",
       x="Chick age (days)", y="Body mass (grams)") +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=13,face="bold"),
        legend.title = element_text(size=13,face="bold"),
        legend.text=element_text(size=13),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major.y = element_line(),
        panel.grid.major.x = element_line())

#plot 3: ckage (days) versus mass (grams) : 2002
ggplot(data2002,aes(x=ckage, y=weight)) +
  geom_point(size = 2,shape = 1, color = "#D9D9D9") +
  theme_classic() +
  labs(title="",
       x="Chick age (days)", y="Body mass (grams)") +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=13,face="bold"),
        legend.title = element_text(size=13,face="bold"),
        legend.text=element_text(size=13),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major.y = element_line(),
        panel.grid.major.x = element_line()) +
  geom_vline(xintercept = 50, color = "black", linetype = "dashed", size = 1)

###############################################
#               Dataset before day 49         #
###############################################

Threshold <- 50
time2002 <- data2002$ckage[data2002$ckage < Threshold]
obs2002  <- data2002$weight[data2002$ckage < Threshold]
data2002cut <- data.frame(time2002,obs2002)
colnames(data2002cut) <- c("ckage","weight")

#plot 4: ckage (days) versus mass (grams) : 2002
ggplot(data2002cut,aes(x=ckage, y=weight)) +
  geom_point(size=1) +
  theme_classic() +
  labs(title="",
       x="Chick age (days)", y="Body mass (grams)") +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=13,face="bold"),
        legend.title = element_text(size=13,face="bold"),
        legend.text=element_text(size=13),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major.y = element_line(),
        panel.grid.major.x = element_line())

###############################################
#             Likelihood approach             #
###############################################
#MLE
Ui_ord<-rbind(c(1,0,0),
              c(0,1,0),
              c(0,0,1))
Ci_ord<-c(1,1,0)
Ui<-rbind(c(1,0,0,0),
          c(0,1,0,0),
          c(0,0,1,0),
          c(0,0,0,1))
Ci<-c(1,1,0,0)
OptMLE_Conf <-constrOptim(c(2,1200,0.5,0.75),FunctionNegativeLogLikFDE_Conformable,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                    outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight)
OptMLE_Cap <-constrOptim(c(2,1200,0.5,0.75),FunctionNegativeLogLikFDE_Caputo,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                    outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,h,m0)
OptMLE_Ord <- constrOptim(c(2,2,0.5),FunctionNegativeLogLikODE,NULL,ui=Ui_ord,ci=Ci_ord,mu=1e-04,control=list(),method="Nelder-Mead",
                          outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight)


saveRDS(OptMLE_Conf,file = "OptMLE_Conf.RDS") 
saveRDS(OptMLE_Cap,file = "OptMLE_Cap.RDS") 
saveRDS(OptMLE_Cap,file = "OptMLE_Cap.RDS") 

###############################################
#             Fitting models                  #
###############################################
timevalues<-seq(min(data2002cut$ckage),max(data2002cut$ckage),length.out=100)

FitModelFDE_Conf<-c()
FitModelODE<-c()
FitModelFDE_Conf_sigma <- c()
FitModelFDE_Caputo<-FunctionPopulationGrowthFDE_Caputo(OptMLE_Cap$par,h,m0)

PopulationGrowthModel_Cap_sigma <- c()
for(i in 1:length(timevalues))
{
  FitModelFDE_Conf[i]<-FunctionPopulationGrowthFDE_Conformable(timevalues[i],OptMLE_Conf$par)
}

for(i in 1:length(timevalues))
{
  FitModelODE[i]<-FunctionPopulationGrowthODE(timevalues[i],OptMLE_Ord$par)
}


for(i in 1:length(data2002cut$ckage))
{
  FitModelFDE_Conf_sigma[i]<-FunctionPopulationGrowthFDE_Conformable(data2002cut$ckage[i],OptMLE_Conf$par)
}

vtau_Conf <- length(data2002cut$ckage)/sum(((data2002cut$weight-FitModelFDE_Conf_sigma)/FitModelFDE_Conf_sigma)^2)
sigma_Conf <- (vtau_Conf)^(-1/2)*FitModelFDE_Conf

for (i in 1:length(data2002cut$ckage)) {
  idx <- which.min(abs(time - data2002cut$ckage[i]))
  PopulationGrowthModel_Cap_sigma[i] <- FitModelFDE_Caputo[idx]
}

vtau_Cap <- length(data2002cut$ckage)/sum(((data2002cut$weight-PopulationGrowthModel_Cap_sigma)/PopulationGrowthModel_Cap_sigma)^2)
sigma_Cap <- (vtau_Cap)^(-1/2)*FitModelFDE_Caputo



Model_ODE  <- data.frame(x = timevalues, y = FitModelODE)
Model_Conf <- data.frame(x = timevalues, y = FitModelFDE_Conf)
Model_Conf_band_1 <- data.frame(x = timevalues, y = FitModelFDE_Conf + 1.96*sigma_Conf)
Model_Conf_band_2 <- data.frame(x = timevalues, y = FitModelFDE_Conf - 1.96*sigma_Conf)
Model_Cap <- data.frame(x = time[141:length(time)], y = FitModelFDE_Caputo[141:length(time)])
Model_Cap_band_1 <- data.frame(x = time[141:length(time)], y = FitModelFDE_Caputo[141:length(time)] + 1.96*sigma_Cap[141:length(time)])
Model_Cap_band_2 <- data.frame(x = time[141:length(time)], y = FitModelFDE_Caputo[141:length(time)] - 1.96*sigma_Cap[141:length(time)])
#plot 4: ckage (days) versus mass (grams) : 2002
ggplot(data2002cut,aes(x=ckage, y=weight)) +
  geom_point(size = 2,shape = 1, color = "#C8C8C8") +
  theme_classic() +
  labs(title="",
       x="Chick age (days)", y="Body mass (grams)") +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=13,face="bold"),
        legend.title = element_text(size=13,face="bold"),
        legend.text=element_text(size=13),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major.y = element_line(),
        panel.grid.major.x = element_line()) +
  geom_line(data = Model_ODE, aes(x = x, y = y), color = "black",linewidth = 1, alpha = 1, linetype = "dotted") + 
  geom_line(data = Model_Conf, aes(x = x, y = y), color = "black",linewidth = 1 ) +
  geom_line(data = Model_Conf_band_1, aes(x = x, y = y), color = "black",linewidth = 1, linetype = "dashed") +
  geom_line(data = Model_Conf_band_2, aes(x = x, y = y), color = "black",linewidth = 1, linetype = "dashed") +
  geom_line(data = Model_Cap, aes(x = x, y = y), color = "gray",linewidth = 1 ) +
  geom_line(data = Model_Cap_band_1, aes(x = x, y = y), color = "gray",linewidth = 1, linetype = "dashed") +
  geom_line(data = Model_Cap_band_2, aes(x = x, y = y), color = "gray",linewidth = 1, linetype = "dashed") 

###############################################
#Likelihood approach: profile alpha           #
###############################################
#restrictions
Ui<-rbind(c(1,0,0),
          c(-1,1,0),
          c(0,0,1))
Ci<-c(0,0,0)

#Profile likelihood alpha: 
valalpha1_Conf<-seq(OptMLE_Conf$par[4],1,length.out=50)
LogProfileLikalpha1_Conf<-c()
vpar0<-OptMLE_Conf$par[-4]
constrOptim(vpar0,FunctionNegativeLogProfileLikFDE_Conformable,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
            outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,valalpha1_Conf[1])
for(i in 1:length(valalpha1_Conf))
{
  OptMLEalpha1<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDE_Conformable,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                            outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,valalpha1_Conf[i])
  LogProfileLikalpha1_Conf[i]<--OptMLEalpha1$value
  vpar0<-OptMLEalpha1$par
  print(i)
}
Relativaalpha1_Conf<-exp(LogProfileLikalpha1_Conf-max(LogProfileLikalpha1_Conf))
plot(valalpha1_Conf,Relativaalpha1_Conf,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(alpha),ylab="Relative profile likelihood",
     main="",ylim=c(0,1))

#Profile likelihood alpha:
valalpha2_Conf<-seq(OptMLE_Conf$par[4],2.9,length.out=50)
LogProfileLikalpha2_Conf<-c()
vpar0<-OptMLE_Conf$par[-4]
for(i in 1:length(valalpha2_Conf))
{
  OptMLEalpha2<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDE_Conformable,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                            outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,valalpha2_Conf[i])
  LogProfileLikalpha2_Conf[i]<--OptMLEalpha2$value
  vpar0<-OptMLEalpha2$par
  print(i)
}
Relativaalpha2_Conf<-exp(LogProfileLikalpha2_Conf-max(LogProfileLikalpha2_Conf))
plot(valalpha2_Conf,Relativaalpha2_Conf,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(alpha),ylab="Relative profile likelihood",
     main="",ylim=c(0,1))


#Plot: Profile alpha
valalpha_Conf<-c(rev(valalpha1_Conf),valalpha2_Conf)
LogProfileLikalpha_Conf<-c(rev(LogProfileLikalpha1_Conf),LogProfileLikalpha2_Conf)
Relativaalpha_Conf<-exp(LogProfileLikalpha_Conf-max(LogProfileLikalpha_Conf))
LimitesIVC_alpha_Conf <- FunctionProfileLCI(Relativaalpha_Conf,valalpha_Conf,0.146)
saveRDS(valalpha_Conf,file = "valalpha_Conf.RDS")
saveRDS(Relativaalpha_Conf,file = "Relativaalpha_Conf.RDS")
plot(valalpha_Conf,Relativaalpha_Conf,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(alpha),ylab="Relative profile likelihood",
     main="",ylim=c(0,1))
segments(OptMLE_Conf$par[4],0,OptMLE_Conf$par[4],1, lty=2)
points(OptMLE_Conf$par[4],0,pch=19,col=1)
text(LimitesIVC_alpha_Conf[1],0,"[",cex=1.5)
text(LimitesIVC_alpha_Conf[2],0,"]",cex=1.5)
L_alpha_Conf <- diff(LimitesIVC_alpha_Conf)

###############################################
#Likelihood approach: profile alpha   Cap     #
###############################################
#restrictions
Ui<-rbind(c(1,0,0),
          c(-1,1,0),
          c(0,0,1))
Ci<-c(0,0,0)

#Profile likelihood alpha: Parte 1 - 30 segundos
valalpha1_Cap<-seq(OptMLE_Cap$par[4],1,length.out=60)
LogProfileLikalpha1_Cap<-c()
vpar0<-c(OptMLE_Cap$par[1],OptMLE_Cap$par[2],OptMLE_Cap$par[3])
constrOptim(vpar0,FunctionNegativeLogProfileLikFDE_Caputo,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
            outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,h,m0,valalpha1_Cap[1])
for(i in 1:length(valalpha1_Cap))
{
  OptMLEalpha1<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDE_Caputo,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                            outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,h,m0,valalpha1_Cap[i])
  LogProfileLikalpha1_Cap[i]<--OptMLEalpha1$value
  vpar0<-c(OptMLE_Cap$par[1],OptMLEalpha1$par[2],OptMLEalpha1$par[3])
  print(i)
}
Relativaalpha1_Cap<-exp(LogProfileLikalpha1_Cap-max(LogProfileLikalpha1_Cap))
plot(valalpha1_Cap,Relativaalpha1_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(alpha),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))

#Profile likelihood alpha: Parte 2 
valalpha2_Cap<-seq(OptMLE_Cap$par[4],1.52,length.out=20)
LogProfileLikalpha2_Cap<-c()
vpar0<-OptMLE_Cap$par[-4]
for(i in 1:length(valalpha2_Cap))
{
  OptMLEalpha2<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDE_Caputo,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                            outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,h,m0,valalpha2_Cap[i])
  LogProfileLikalpha2_Cap[i]<--OptMLEalpha2$value
  vpar0<-c(OptMLE_Cap$par[1],OptMLEalpha2$par[2],OptMLEalpha2$par[3])
  print(-OptMLEalpha2$value)
  print(i)
}
Relativaalpha2_Cap<-exp(LogProfileLikalpha2_Cap-max(LogProfileLikalpha2_Cap))
plot(valalpha2_Cap,Relativaalpha2_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(alpha),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))

#Profile likelihood alpha: Parte 3 
valalpha3_Cap<-seq(1.65,1.52,length.out=10)
LogProfileLikalpha3_Cap<-c()
vpar0<-c(15,OptMLEalpha2$par[2],OptMLEalpha2$par[3])
for(i in 1:length(valalpha3_Cap))
{
  OptMLEalpha3<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDE_Caputo,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                            outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,h,m0,valalpha3_Cap[i])
  LogProfileLikalpha3_Cap[i]<--OptMLEalpha3$value
  #vpar0<-c(10,OptMLEalpha2$par[2],OptMLEalpha2$par[3])
  print(-OptMLEalpha3$value)
  print(i)
}
Relativaalpha3_Cap<-exp(LogProfileLikalpha3_Cap-max(LogProfileLikalpha3_Cap))
plot(valalpha3_Cap,Relativaalpha3_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(alpha),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))



#Plot: Profile alpha
valalpha_Cap<-c(rev(valalpha1_Cap),valalpha2_Cap,rev(valalpha3_Cap))
LogProfileLikalpha_Cap<-c(rev(LogProfileLikalpha1_Cap),LogProfileLikalpha2_Cap,rev(LogProfileLikalpha3_Cap))
Relativaalpha_Cap<-exp(LogProfileLikalpha_Cap-max(LogProfileLikalpha_Cap))
LimitesIVC_alpha_Cap <- FunctionProfileLCI(Relativaalpha_Cap,valalpha_Cap,0.146)
saveRDS(valalpha_Cap,file = "valalpha_Cap.RDS")
saveRDS(Relativaalpha_Cap,file = "Relativaalpha_Cap.RDS")
plot(valalpha_Cap,Relativaalpha_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(alpha),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))
segments(OptMLE_Cap$par[4],0,OptMLE_Cap$par[4],1, lty=2)
points(OptMLE_Cap$par[4],0,pch=19,col=1)
text(LimitesIVC_alpha_Cap[1],0,"[",cex=1.5)
text(LimitesIVC_alpha_Cap[2],0,"]",cex=1.5)
L_alpha_Cap <- diff(LimitesIVC_alpha_Cap)

#Both
df_alpha <- rbind(
  data.frame(x = valalpha_Cap, y = Relativaalpha_Cap, curva = "Caputo"),
  data.frame(x = valalpha_Conf, y = Relativaalpha_Conf, curva = "Conformable")
)
# Graficar
ggplot(df_alpha, aes(x = x, y = y, color = curva)) +
  geom_line(size = 1) +
  geom_segment(aes(x = LimitesIVC_alpha_Cap[1], xend = LimitesIVC_alpha_Cap[2], 
                   y = 0.146, yend = 0.146),
               linetype = "dashed", color = "gray", size = 1) +
  geom_segment(aes(x = LimitesIVC_alpha_Conf[1], xend = LimitesIVC_alpha_Conf[2], 
                   y = 0.146, yend = 0.146),
               linetype = "dashed", color = "black", size = 1) +
  labs(title = "", x = expression(alpha), y = "Relative profile likelihood") +
  scale_color_manual(values = c("gray", "black")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) 

###############################################
#Likelihood approach: profile x0              #
###############################################
#restrictions
Ui<-rbind(c(1,0,0),
          c(0,1,0),
          c(0,0,1))
#Profile likelihood x0:
valx0_Conf<-seq(150,400,length.out=50)
LogProfileLikx0_Conf<-c()
vpar0<-OptMLE_Conf$par[-1]
for(i in 1:length(valx0_Conf))
{
  Ci<-c(valx0_Conf[i],0,0)
  OptMLEx0<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDEx0_Conformable,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                            outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,valx0_Conf[i])
  LogProfileLikx0_Conf[i]<--OptMLEx0$value
  vpar0<-OptMLEx0$par
  print(i)
}
Relativax0_Conf<-exp(LogProfileLikx0_Conf-max(LogProfileLikx0_Conf))
saveRDS(valx0_Conf,file = "valx0_Conf.RDS")
saveRDS(Relativax0_Conf,file = "Relativax0_Conf.RDS")
LimitesIVC_x0_Conf <- FunctionProfileLCI(Relativax0_Conf,valx0_Conf,0.146)
plot(valx0_Conf,Relativax0_Conf,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(x[0]),ylab="Relative profile likelihood",
     main="",ylim=c(0,1))
segments(OptMLE_Conf$par[1],0,OptMLE_Conf$par[1],1, lty=2)
points(OptMLE_Conf$par[1],0,pch=19,col=1)
text(LimitesIVC_x0_Conf[1],0,"[",cex=1.5)
text(LimitesIVC_x0_Conf[2],0,"]",cex=1.5)
L_x0_Conf <- diff(LimitesIVC_x0_Conf)


###############################################
#Likelihood approach: profile x0    Caputo    #
###############################################
#restrictions
Ui<-rbind(c(1,0,0),
          c(0,1,0),
          c(0,0,1))

#Profile likelihood x0: Parte 1
valx01_Cap<-seq(OptMLE_Cap$par[1],2,length.out=50) 
LogProfileLikx01_Cap<-c()
vpar0<-OptMLE_Cap$par[-1]
for(i in 1:length(valx01_Cap))
{
  Ci<-c(valx01_Cap[i],0,0)
  OptMLEx0<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDEx0_Caputo,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                        outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,h,m0,valx01_Cap[i])
  LogProfileLikx01_Cap[i]<--OptMLEx0$value
  vpar0<-OptMLEx0$par
  print(i)
}
Relativax01_Cap<-exp(LogProfileLikx01_Cap-max(LogProfileLikx01_Cap))
plot(valx01_Cap,Relativax01_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(x[0]),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))

#Profile likelihood x0: Parte 2
valx02_Cap<-seq(OptMLE_Cap$par[1],80,length.out=50) 
LogProfileLikx02_Cap<-c()
vpar0<-OptMLE_Cap$par[-1]
for(i in 1:length(valx02_Cap))
{
  Ci<-c(valx02_Cap[i],0,0)
  OptMLEx0<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDEx0_Caputo,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                        outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,h,m0,valx02_Cap[i])
  LogProfileLikx02_Cap[i]<--OptMLEx0$value
  vpar0<-OptMLEx0$par
  print(i)
}
Relativax02_Cap<-exp(LogProfileLikx02_Cap-max(LogProfileLikx02_Cap))
plot(valx02_Cap,Relativax02_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(x[0]),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))

#Plot: Profile x0
valx0_Cap<-c(rev(valx01_Cap),valx02_Cap)
LogProfileLikx0_Cap<-c(rev(LogProfileLikx01_Cap),LogProfileLikx02_Cap)
Relativax0_Cap<-exp(LogProfileLikx0_Cap-max(LogProfileLikx0_Cap))
LimitesIVC_x0_Cap <- FunctionProfileLCI(Relativax0_Cap,valx0_Cap,0.146)
saveRDS(valx0_Cap,file = "valx0_Cap.RDS")
saveRDS(Relativax0_Cap,file = "Relativax0_Cap.RDS")
plot(valx0_Cap,Relativax0_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(x[0]),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))
segments(OptMLE_Cap$par[1],0,OptMLE_Cap$par[1],1, lty=2)
points(OptMLE_Cap$par[1],0,pch=19,col=1)
text(LimitesIVC_x0_Cap[1],0,"[",cex=1.5)
text(LimitesIVC_x0_Cap[2],0,"]",cex=1.5)
L_x0_Cap <- diff(LimitesIVC_x0_Cap)

#Both
df_x0 <- rbind(
  data.frame(x = valx0_Cap, y = Relativax0_Cap, curva = "Caputo"),
  data.frame(x = valx0_Conf, y = Relativax0_Conf, curva = "Conformable"))

ggplot(df_x0, aes(x = x, y = y, color = curva)) +
  geom_line(size = 1) +
  geom_segment(aes(x = LimitesIVC_x0_Cap[1], xend = LimitesIVC_x0_Cap[2], 
                   y = 0.146, yend = 0.146),
               linetype = "dashed", color = "gray", size = 1) +
  geom_segment(aes(x = LimitesIVC_x0_Conf[1], xend = LimitesIVC_x0_Conf[2], 
                   y = 0.146, yend = 0.146),
               linetype = "dashed", color = "black", size = 1) +
  labs(title = "", x = expression(x[0]), y = "Relative profile likelihood") +
  scale_color_manual(values = c("gray", "black")) +
  theme_minimal() +
  theme(legend.position = "none") 



###############################################
#Likelihood approach: profile K               #
###############################################
#restrictions
Ui<-rbind(c(-1,0,0),
          c(0,1,0),
          c(0,0,1))
#Profile likelihood x0:
valK_Conf<-seq(1040,1090,length.out=50)
LogProfileLikK_Conf<-c()
vpar0<-OptMLE_Conf$par[-2]
for(i in 1:length(valK_Conf))
{
  Ci<-c(-valK_Conf[i],0,0)
  OptMLEK<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDEK_Conformable,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                        outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,valK_Conf[i])
  LogProfileLikK_Conf[i]<--OptMLEK$value
  vpar0<-OptMLEK$par
  print(i)
}
RelativaK_Conf<-exp(LogProfileLikK_Conf-max(LogProfileLikK_Conf))
saveRDS(valK_Conf,file = "valK_Conf.RDS")
saveRDS(RelativaK_Conf,file = "RelativaK_Conf.RDS")
LimitesIVC_K_Conf <- FunctionProfileLCI(RelativaK_Conf,valK_Conf,0.146)
plot(valK_Conf,RelativaK_Conf,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(K),ylab="Relative profile likelihood",
     main="",ylim=c(0,1))
segments(OptMLE_Conf$par[2],0,OptMLE_Conf$par[2],1, lty=2)
points(OptMLE_Conf$par[2],0,pch=19,col=1)
text(LimitesIVC_K_Conf[1],0,"[",cex=1.5)
text(LimitesIVC_K_Conf[2],0,"]",cex=1.5)
L_K <- diff(LimitesIVC_K_Conf)


###############################################
#Likelihood approach: profile K   Caputo      #
###############################################
#restrictions
Ui<-rbind(c(-1,0,0),
          c(0,1,0),
          c(0,0,1))
#Profile likelihood x0:
valK_Cap<-seq(750,1060,length.out=50)
LogProfileLikK_Cap<-c()
vpar0<-OptMLE_Cap$par[-2]
for(i in 1:length(valK_Cap))
{
  Ci<-c(-valK_Cap[i],0,0)
  OptMLEK<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDEK_Caputo,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                       outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,h,m0,valK_Cap[i])
  LogProfileLikK_Cap[i]<--OptMLEK$value
  vpar0<-OptMLEK$par
  print(i)
}
RelativaK_Cap<-exp(LogProfileLikK_Cap-max(LogProfileLikK_Cap))
saveRDS(valK_Cap,file = "valK_Cap.RDS")
saveRDS(RelativaK_Cap,file = "RelativaK_Cap.RDS")
LimitesIVC_K_Cap <- FunctionProfileLCI(RelativaK_Cap,valK_Cap,0.146)
plot(valK_Cap,RelativaK_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(K),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))
segments(OptMLE_Cap$par[2],0,OptMLE_Cap$par[2],1, lty=2)
points(OptMLE_Cap$par[2],0,pch=19,col=1)
text(LimitesIVC_K_Cap[1],0,"[",cex=1.5)
text(LimitesIVC_K_Cap[2],0,"]",cex=1.5)
L_K_Cap <- diff(LimitesIVC_K_Cap)

#Both
df_K <- rbind(
  data.frame(x = valK_Cap, y = RelativaK_Cap, curva = "Caputo"),
  data.frame(x = valK_Conf, y = RelativaK_Conf, curva = "Conformable"))

ggplot(df_K, aes(x = x, y = y, color = curva)) +
  geom_line(size = 1) +
  geom_segment(aes(x = LimitesIVC_K_Cap[1], xend = LimitesIVC_K_Cap[2], 
                   y = 0.146, yend = 0.146),
               linetype = "dashed", color = "gray", size = 1) +
  geom_segment(aes(x = LimitesIVC_K_Conf[1], xend = LimitesIVC_K_Conf[2], 
                   y = 0.146, yend = 0.146),
               linetype = "dashed", color = "black", size = 1) +
  labs(title = "", x = expression(K), y = "Relative profile likelihood") +
  scale_color_manual(values = c("gray", "black")) +
  theme_minimal() +
  theme(legend.position = "none") 

###############################################
#Likelihood approach: profile lambda          #
###############################################
#restrictions
Ui<-rbind(c(1,0,0),
          c(-1,1,0),
          c(0,0,1))
Ci<-c(0,0,0)
#Profile likelihood lambda
vallambda1_Conf<-seq(OptMLE_Conf$par[3],0.0001,length.out=50)
LogProfileLiklambda1_Conf<-c()
vpar0<-OptMLE_Conf$par[-3]
for(i in 1:length(vallambda1_Conf))
{
  OptMLElambda1<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDElambda_Conformable,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                       outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,vallambda1_Conf[i])
  LogProfileLiklambda1_Conf[i]<--OptMLElambda1$value
  vpar0<-OptMLElambda1$par
  print(i)
}
Relativalambda1_Conf<-exp(LogProfileLiklambda1_Conf-max(LogProfileLiklambda1_Conf))
plot(vallambda1_Conf,Relativalambda1_Conf,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(lambda),ylab="Relative profile likelihood",
     main="",ylim=c(0,1))
segments(OptMLE_Conf$par[3],0,OptMLE_Conf$par[3],1, lty=2)
points(OptMLE_Conf$par[3],0,pch=19,col=1)

#Profile likelihood lambda
vallambda2_Conf<-seq(OptMLE_Conf$par[3],0.07,length.out=50)
LogProfileLiklambda2_Conf<-c()
vpar0<-OptMLE_Conf$par[-3]
for(i in 1:length(vallambda2_Conf))
{
  OptMLElambda2<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDElambda_Conformable,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                             outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,vallambda2_Conf[i])
  LogProfileLiklambda2_Conf[i]<--OptMLElambda2$value
  vpar0<-OptMLElambda2$par
  print(i)
}
Relativalambda2_Conf<-exp(LogProfileLiklambda2_Conf-max(LogProfileLiklambda2_Conf))
plot(vallambda2_Conf,Relativalambda2_Conf,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(lambda),ylab="Relative profile likelihood",
     main="",ylim=c(0,1))
segments(OptMLE_Conf$par[3],0,OptMLE_Conf$par[3],1, lty=2)
points(OptMLE_Conf$par[3],0,pch=19,col=1)

#Plot: Profile lambda
vallambda_Conf<-c(rev(vallambda1_Conf),vallambda2_Conf)
LogProfileLiklambda_Conf<-c(rev(LogProfileLiklambda1_Conf),LogProfileLiklambda2_Conf)
Relativalambda_Conf<-exp(LogProfileLiklambda_Conf-max(LogProfileLiklambda_Conf))
saveRDS(vallambda_Conf,file = "vallambda_Conf.RDS")
saveRDS(Relativalambda_Conf,file = "Relativalambda_Conf.RDS")
LimitesIVC_lambda_Conf <- FunctionProfileLCI(Relativalambda_Conf,vallambda_Conf,0.146)
plot(vallambda_Conf,Relativalambda_Conf,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(lambda),ylab="Relative profile likelihood",
     main="",ylim=c(0,1))
segments(OptMLE_Conf$par[3],0,OptMLE_Conf$par[3],1, lty=2)
points(OptMLE_Conf$par[3],0,pch=19,col=1)
text(LimitesIVC_lambda_Conf[1],0,"[",cex=1.5)
text(LimitesIVC_lambda_Conf[2],0,"]",cex=1.5)
L_lambda <- diff(LimitesIVC_lambda_Conf)



###############################################
#Likelihood approach: profile lambda   Caputo #
###############################################
#restrictions
Ui<-rbind(c(1,0,0),
          c(-1,1,0),
          c(0,0,1))
Ci<-c(0,0,0)
#Profile likelihood lambda: 
vallambda1_Cap<-seq(OptMLE_Cap$par[3],0.005,length.out=50)
LogProfileLiklambda1_Cap<-c()
vpar0<-OptMLE_Cap$par[-3]
for(i in 1:length(vallambda1_Cap))
{
  OptMLElambda1<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDElambda_Caputo,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                             outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,h,m0,vallambda1_Cap[i])
  LogProfileLiklambda1_Cap[i]<--OptMLElambda1$value
  vpar0<-OptMLElambda1$par
  print(i)
}
Relativalambda1_Cap<-exp(LogProfileLiklambda1_Cap-max(LogProfileLiklambda1_Cap))
plot(vallambda1_Cap,Relativalambda1_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(lambda),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))
segments(OptMLE_Cap$par[3],0,OptMLE_Cap$par[3],1, lty=2)
points(OptMLE_Cap$par[3],0,pch=19,col=1)

#Profile likelihood lambda
vallambda2_Cap<-seq(OptMLE_Cap$par[3],0.1,length.out=50)
LogProfileLiklambda2_Cap<-c()
vpar0<-OptMLE_Cap$par[-3]
for(i in 1:length(vallambda2_Cap))
{
  OptMLElambda2<-constrOptim(vpar0,FunctionNegativeLogProfileLikFDElambda_Caputo,NULL,ui=Ui,ci=Ci,mu=1e-04,control=list(),method="Nelder-Mead",
                             outer.iterations = 100, outer.eps = 1e-05,data2002cut$ckage,data2002cut$weight,h,m0,vallambda2_Cap[i])
  LogProfileLiklambda2_Cap[i]<--OptMLElambda2$value
  vpar0<-OptMLElambda2$par
  print(i)
}
Relativalambda2_Cap<-exp(LogProfileLiklambda2_Cap-max(LogProfileLiklambda2_Cap))
plot(vallambda2_Cap,Relativalambda2_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(lambda),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))
segments(OptMLE_Cap$par[3],0,OptMLE_Cap$par[3],1, lty=2)
points(OptMLE_Cap$par[3],0,pch=19,col=1)

#Plot: Profile lambda
vallambda_Cap<-c(rev(vallambda1_Cap),vallambda2_Cap)
LogProfileLiklambda_Cap<-c(rev(LogProfileLiklambda1_Cap),LogProfileLiklambda2_Cap)
Relativalambda_Cap<-exp(LogProfileLiklambda_Cap-max(LogProfileLiklambda_Cap))
saveRDS(vallambda_Cap,file = "vallambda_Cap.RDS")
saveRDS(Relativalambda_Cap,file = "Relativalambda_Cap.RDS")
LimitesIVC_lambda_Cap <- FunctionProfileLCI(Relativalambda_Cap,vallambda_Cap,0.146)
plot(vallambda_Cap,Relativalambda_Cap,type="l",lwd=2,col=1,cex.lab=1.20,
     xlab=expression(lambda),ylab="Relative profile likelihood Caputo",
     main="",ylim=c(0,1))
segments(OptMLE_Cap$par[3],0,OptMLE_Cap$par[3],1, lty=2)
points(OptMLE_Cap$par[3],0,pch=19,col=1)
text(LimitesIVC_lambda_Cap[1],0,"[",cex=1.5)
text(LimitesIVC_lambda_Cap[2],0,"]",cex=1.5)
L_lambda_Cap <- diff(LimitesIVC_lambda_Cap)

#Both
df_lambda <- rbind(
  data.frame(x = vallambda_Cap, y = Relativalambda_Cap, curva = "Caputo"),
  data.frame(x = vallambda_Conf, y = Relativalambda_Conf, curva = "Conformable"))

ggplot(df_lambda, aes(x = x, y = y, color = curva)) +
  geom_line(size = 1) +
  geom_segment(aes(x = LimitesIVC_lambda_Cap[1], xend = LimitesIVC_lambda_Cap[2], 
                   y = 0.146, yend = 0.146),
               linetype = "dashed", color = "gray", size = 1) +
  geom_segment(aes(x = LimitesIVC_lambda_Conf[1], xend = LimitesIVC_lambda_Conf[2], 
                   y = 0.146, yend = 0.146),
               linetype = "dashed", color = "black", size = 1) +
  labs(title = "", x = expression(lambda), y = "Relative profile likelihood") +
  scale_color_manual(values = c("gray", "black")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) 


###############################################
#               Efficiency                    #
###############################################
Estimated_Ord <- c()
n<-length(data2002cut$ckage)
for(i in 1:n)
{
  Estimated_Ord[i]<- FunctionPopulationGrowthODE(data2002cut$ckage[i],OptMLE_Ord$par)
}
E_Ord <- sum(((data2002cut$weight-Estimated_Ord))^2)

Estimated_Conf <- c()
n<-length(data2002cut$ckage)
for(i in 1:n)
{
  Estimated_Conf[i]<- FunctionPopulationGrowthFDE_Conformable(data2002cut$ckage[i],OptMLE_Conf$par)
}
E_Conf <- sum(((data2002cut$weight-Estimated_Conf))^2)

Estimated_Cap <- c()
Sol_Cap <- FunctionPopulationGrowthFDE_Caputo(OptMLE_Cap$par,h,m0)
for (i in 1:n) {
  idx <- which.min(abs(time - data2002cut$ckage[i]))
  Estimated_Cap[i] <- Sol_Cap[idx]
}
E_Cap <- sum(((data2002cut$weight-Estimated_Cap))^2)


Eff_Conf <- abs((E_Ord-E_Conf)/E_Ord)
Eff_Cap <- abs((E_Ord-E_Cap)/E_Ord)






###############################################
#               Prediction                    #
###############################################
timevalues<-seq(min(data2002$ckage),max(data2002$ckage),length.out=100)
m0 <- 720
final_time <- h*m0
time_2 <- seq(0,final_time,h)
FitModelFDE_Conf<-c()
FitModelODE<-c()
FitModelFDE_Conf_sigma <- c()
FitModelFDE_Caputo<-FunctionPopulationGrowthFDE_Caputo(OptMLE_Cap$par,h,m0)

PopulationGrowthModel_Cap_sigma <- c()
for(i in 1:length(timevalues))
{
  FitModelFDE_Conf[i]<-FunctionPopulationGrowthFDE_Conformable(timevalues[i],OptMLE_Conf$par)
}

for(i in 1:length(timevalues))
{
  FitModelODE[i]<-FunctionPopulationGrowthODE(timevalues[i],OptMLE_Ord$par)
}


for(i in 1:length(data2002$ckage))
{
  FitModelFDE_Conf_sigma[i]<-FunctionPopulationGrowthFDE_Conformable(data2002$ckage[i],OptMLE_Conf$par)
}

vtau_Conf <- length(data2002$ckage)/sum(((data2002$weight-FitModelFDE_Conf_sigma)/FitModelFDE_Conf_sigma)^2)
sigma_Conf <- (vtau_Conf)^(-1/2)*FitModelFDE_Conf

for (i in 1:length(data2002$ckage)) {
  idx <- which.min(abs(time_2 - data2002$ckage[i]))
  PopulationGrowthModel_Cap_sigma[i] <- FitModelFDE_Caputo[idx]
}

vtau_Cap <- length(data2002$ckage)/sum(((data2002$weight-PopulationGrowthModel_Cap_sigma)/PopulationGrowthModel_Cap_sigma)^2)
sigma_Cap <- (vtau_Cap)^(-1/2)*FitModelFDE_Caputo



Model_ODE  <- data.frame(x = timevalues, y = FitModelODE)
Model_Conf <- data.frame(x = timevalues, y = FitModelFDE_Conf)
Model_Conf_band_1 <- data.frame(x = timevalues, y = FitModelFDE_Conf + 1.96*sigma_Conf)
Model_Conf_band_2 <- data.frame(x = timevalues, y = FitModelFDE_Conf - 1.96*sigma_Conf)
Model_Cap <- data.frame(x = time_2[141:length(time_2)], y = FitModelFDE_Caputo[141:length(time_2)])
Model_Cap_band_1 <- data.frame(x = time_2[141:length(time_2)], y = FitModelFDE_Caputo[141:length(time_2)] + 1.96*sigma_Cap[141:length(time_2)])
Model_Cap_band_2 <- data.frame(x = time_2[141:length(time_2)], y = FitModelFDE_Caputo[141:length(time_2)] - 1.96*sigma_Cap[141:length(time_2)])
#plot 4: ckage (days) versus mass (grams) : 2002
ggplot(data2002,aes(x=ckage, y=weight)) +
  geom_point(size = 2,shape = 1, color = "#C8C8C8") +
  theme_classic() +
  labs(title="",
       x="Chick age (days)", y="Body mass (grams)") +
  theme(axis.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title=element_text(size=13,face="bold"),
        legend.title = element_text(size=13,face="bold"),
        legend.text=element_text(size=13),
        plot.title = element_text(size = 14, face = "bold"),
        panel.grid.major.y = element_line(),
        panel.grid.major.x = element_line()) +
  geom_line(data = Model_Conf, aes(x = x, y = y), color = "black",linewidth = 1 ) +
  geom_line(data = Model_Conf_band_1, aes(x = x, y = y), color = "black",linewidth = 1, linetype = "dashed") +
  geom_line(data = Model_Conf_band_2, aes(x = x, y = y), color = "black",linewidth = 1, linetype = "dashed") +
  geom_line(data = Model_Cap, aes(x = x, y = y), color = "gray",linewidth = 1 ) +
  geom_line(data = Model_Cap_band_1, aes(x = x, y = y), color = "gray",linewidth = 1, linetype = "dashed") +
  geom_line(data = Model_Cap_band_2, aes(x = x, y = y), color = "gray",linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = Threshold, color = "black", linetype = "dotted", size = 1)
