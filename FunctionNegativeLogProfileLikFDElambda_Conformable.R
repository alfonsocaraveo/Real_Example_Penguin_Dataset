FunctionNegativeLogProfileLikFDElambda_Conformable<-function(VecPar124,VecTimeObs,VecObs,VecPar3)
{
  VecPar<-c(VecPar124[-3],VecPar3,VecPar124[3])
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

