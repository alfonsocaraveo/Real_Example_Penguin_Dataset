FunctionNegativeLogLikFDE_Caputo<-function(VecPar,VecTimeObs,VecObs,h,m0)
{
  PopulationGrowthModel<-c()
  n<-length(VecObs)
  Sol <- FunctionPopulationGrowthFDE_Caputo(VecPar,h,m0)
  final_time<- h*m0
  time <- seq(0,final_time,h)
  for (i in 1:n) {
    idx <- which.min(abs(time - VecTimeObs[i]))
    #print(idx)
    PopulationGrowthModel[i] <- Sol[idx]
    #print(PopulationGrowthModel[i])
  }
  vtau<-n/sum(((VecObs-PopulationGrowthModel)/PopulationGrowthModel)^2)
  LogLik<-(n/2)*log(vtau)-sum(log(PopulationGrowthModel))
  return(-LogLik)
}

