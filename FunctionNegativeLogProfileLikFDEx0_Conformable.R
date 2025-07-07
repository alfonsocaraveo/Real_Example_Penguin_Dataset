FunctionNegativeLogProfileLikFDEx0_Conformable<-function(VecPar234,VecTimeObs,VecObs,VecPar1)
{
  VecPar<-c(VecPar1,VecPar234)
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

