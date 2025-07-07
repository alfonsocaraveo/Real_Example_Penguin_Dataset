FunctionPopulationGrowthODE<-function(TimeValue,VecPar)
{
  vx0<-VecPar[1]
  vK<-VecPar[2]
  vlambda<-VecPar[3]
  vf<-(vK*vx0)/(vx0+(vK-vx0)*exp(-vlambda*(TimeValue)))
  return(vf)
}