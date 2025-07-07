FunctionNegativeLogProfileLikFDEK_Conformable<-function(VecPar134,VecTimeObs,VecObs,VecPar2)
{
  VecPar<-c(VecPar134[1],VecPar2,VecPar134[-1])
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

