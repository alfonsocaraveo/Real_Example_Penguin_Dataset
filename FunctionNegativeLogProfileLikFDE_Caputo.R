FunctionNegativeLogProfileLikFDE_Caputo<-function(VecPar123,VecTimeObs,VecObs,h,m0,VecPar4)
{
  VecPar<-c(VecPar123,VecPar4)
  PopulationGrowthModel<-c()
  n<-length(VecTimeObs)
  Sol <- FunctionPopulationGrowthFDE_Caputo(VecPar,h,m0)
  final_time<- h*m0
  time <- seq(0,final_time,h)
  for (i in 1:n) {
    idx <- which.min(abs(time - VecTimeObs[i]))
    PopulationGrowthModel[i] <- Sol[idx]
    #print(idx)
  }
  vtau<-n/sum(((VecObs-PopulationGrowthModel)/PopulationGrowthModel)^2)
  LogLik<-(n/2)*log(vtau)-sum(log(PopulationGrowthModel))
  return(-LogLik)
}

