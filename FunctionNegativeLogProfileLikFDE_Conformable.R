FunctionNegativeLogProfileLikFDE_Conformable<-function(VecPar123,VecTimeObs,VecObs,VecPar4)
{
  VecPar<-c(VecPar123,VecPar4)
  PopulationGrowthModel<-c()
  n<-length(VecTimeObs)
  for(i in 1:n)
  {
    PopulationGrowthModel[i]<-FunctionPopulationGrowthFDE_Conformable(VecTimeObs[i],VecPar)
  }
  vtau<-n/sum(((VecObs-PopulationGrowthModel)/PopulationGrowthModel)^2)
  LogLik<-(n/2)*log(vtau)-sum(log(PopulationGrowthModel))
  return(-LogLik)
}

