FunctionPopulationGrowthFDE_Caputo <- function(VecPar,h,m0)
{
  omega_k <- c()
  omega_k[1] <- (1-(1+VecPar[4])/(1))*1
  for (j in 1:m0) {
    omega_k[j+1] <- (1-(1+VecPar[4])/(j+1))*omega_k[j]
  }
  
  
  Solucion <- c()
  Solucion[1] <- VecPar[1]
  
  for (m in 1:m0) {
    
    
    suma <- 0
    for (k in 1:m) {
      suma <- suma + omega_k[k]*Solucion[m-k+1]
    }
    
    Solucion[m+1] <- h^(VecPar[4])*VecPar[3]*Solucion[m]*(1-Solucion[m]/VecPar[2])-suma
    
    
    suma <- 0
  }
  return(Solucion)
}