FunctionNegativeLogLikODE<-function(VecPar,VecTimeObs,VecObs)
{
  PopulationGrowthModel<-c()
  n<-length(VecTimeObs)
  for(i in 1:n)
  {
    PopulationGrowthModel[i]<-FunctionPopulationGrowthODE(VecTimeObs[i],VecPar)
  }
  vtau<-n/sum(((VecObs-PopulationGrowthModel)/PopulationGrowthModel)^2)
  LogLik<-(n/2)*log(vtau)-sum(log(PopulationGrowthModel))
  return(-LogLik)
}

