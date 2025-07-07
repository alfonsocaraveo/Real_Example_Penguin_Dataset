FunctionNegativeLogProfileLikFDEK_Caputo<-function(VecPar134,VecTimeObs,VecObs,h,m0,VecPar2)
{
  VecPar<-c(VecPar134[1],VecPar2,VecPar134[-1])
  PopulationGrowthModel<-c()
  n<-length(VecTimeObs)
  Sol <- FunctionPopulationGrowthFDE_Caputo(VecPar,h,m0)
  final_time<- h*m0
  time <- seq(0,final_time,h)
  for (i in 1:n) {
    idx <- which.min(abs(time - VecTimeObs[i]))
    PopulationGrowthModel[i] <- Sol[idx]
    #print(PopulationGrowthModel[i])
  }
  vtau<-n/sum(((VecObs-PopulationGrowthModel)/PopulationGrowthModel)^2)
  LogLik<-(n/2)*log(vtau)-sum(log(PopulationGrowthModel))
  return(-LogLik)
}
