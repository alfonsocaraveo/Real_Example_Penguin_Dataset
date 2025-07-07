# Repository containing model selection procedures applied to a penguin dataset
Repository containing model selection procedures applied to the penguin dataset, 
based on the article by A. Chiaradia and I. C. Nisbet on plasticity in parental provisioning and 
chick growth in little penguins (Eudyptula minor) during years of high and low breeding success.

# Scripts Overview for the model selection analysis
This repository contains the following R scripts, each serving a specific purpose in the workflow:

# 1. .R
**Description:** Script for model selection, where we estimate the MLE of the Caputo, conformable, and classical logistic models using data up to day 49.
We plot chick age versus body mass growth for the year 2002, and provide a graphical comparison of the models, along with the calculation of normal bands for each one.

Additionally, we compute the relative profile likelihood for each parameter **alpha**, **x0**, **K**, and **lambda**.
For each parameter, we calculate the likelihood-based confidence intervals and perform both graphical and numerical comparisons.

We also compute the efficiency index for the Caputo and conformable models. Finally, we assess the predictive performance of each model beyond day 49. 

**Libraries Used:** 

- `nlme`: Used for example datasets.
- `FlexParamCurve`: Used for example datasets.

- `ggplot2`: Used for plotting.
- `ggbreak`: Used for plotting.
- `ggthemes`: Used for plotting.
- `fields`: Used for plotting.

# 2. FunctionPopulationGrowthFDE_Conformable.R
**Description:** Implements the solution of the fractional logistic model considering the conformable derivative.

**Input:**
- `TimeValue`: A vector of time points associated with each observation.
- `VecPar`: A vector containing the parameters **x0**, **K**, **lambda**, and **alpha**.


# 3. FunctionPopulationGrowthODE.R
**Description:** Implements the solution of the classical logistic model.

**Input:**
- `TimeValue`: A vector of time points associated with each observation.
- `VecPar`: A vector containing the parameters **x0**, **K**, and **lambda**.


# 4. FunctionPopulationGrowthFDE_Caputo.R
**Description:** Implements the solution of the fractional logistic model considering the Caputo derivative.

**Input:**
- `VecPar`: A vector containing the parameters **x0**, **K**, and **lambda**.
- `h`: The step size for the computation.
- `m0`: The total number of points to generate.

# 5. FunctionNegativeLogLikFDE_Conformable.R
**Description:** Implements the negative loglikelihood function for the parameters **x0**, **K**, **lambda**, and **alpha** considering the conformable derivative.

**Input:**
- `VecPar`: A vector containing the parameters **x0**, **K**, **lambda**, and **alpha**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.

# 6. FunctionNegativeLogLikFDE_Caputo.R
**Description:** Implements the negative loglikelihood function for the parameters **x0**, **K**, **lambda**, and **alpha** considering the Caputo derivative.

**Input:**
- `VecPar`: A vector containing the parameters **x0**, **K**, **lambda**, and **alpha**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.
- `h`: The step size for the computation.
- `m0`: The total number of points to generate.

# 7. FunctionNegativeLogLikODE.R
**Description:** Implements the negative loglikelihood function for the parameters **x0**, **K**, and **lambda** considering the classical derivative.

**Input:**
- `VecPar`: A vector containing the parameters **x0**, **K**, and **lambda**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.

# 8. FunctionNegativeLogProfileLikFDE_Conformable.R
**Description:** Implements the negative profile loglikelihood function for the the parameter of interest **alpha** considering the conformable derivative.

**Input:**
- `VecPar123`: A vector containing the parameters **x0**, **K**, and **lambda**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.
- `VecPar4`: A numeric value for the parameter of interest, **alpha**.

 # 9. FunctionNegativeLogProfileLikFDE_Caputo.R
**Description:** Implements the negative profile loglikelihood function for the the parameter of interest **alpha** considering the Caputo derivative.

**Input:**
- `VecPar123`: A vector containing the parameters **x0**, **K**, and **lambda**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.
- `h`: The step size for the computation.
- `m0`: The total number of points to generate.
- `VecPar4`: A numeric value for the parameter of interest, **alpha**.

# 10. FunctionNegativeLogProfileLikFDEx0_Conformable.R
**Description:** Implements the negative profile loglikelihood function for the the parameter of interest **x0** considering the conformable derivative.

**Input:**
- `VecPar234`: A vector containing the parameters **K**, **lambda**, and **alpha**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.
- `VecPar1`: A numeric value for the parameter of interest, **x0**.

 # 11. FunctionNegativeLogProfileLikFDEx0_Caputo.R
**Description:** Implements the negative profile loglikelihood function for the the parameter of interest **x0** considering the Caputo derivative.

**Input:**
- `VecPar234`: A vector containing the parameters **K**, **lambda**, and **alpha**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.
- `h`: The step size for the computation.
- `m0`: The total number of points to generate.
- `VecPar1`: A numeric value for the parameter of interest, **x0**.

# 12. FunctionNegativeLogProfileLikFDEK_Conformable.R
**Description:** Implements the negative profile loglikelihood function for the the parameter of interest **K** considering the conformable derivative.

**Input:**
- `VecPar134`: A vector containing the parameters **x0**, **lambda**, and **alpha**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.
- `VecPar2`: A numeric value for the parameter of interest, **K**.

 # 13. FunctionNegativeLogProfileLikFDEK_Caputo.R
**Description:** Implements the negative profile loglikelihood function for the the parameter of interest **K** considering the Caputo derivative.

**Input:**
- `VecPar134`: A vector containing the parameters **x0**, **lambda**, and **alpha**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.
- `h`: The step size for the computation.
- `m0`: The total number of points to generate.
- `VecPar2`: A numeric value for the parameter of interest, **K**.

# 14. FunctionNegativeLogProfileLikFDElambda_Conformable.R
**Description:** Implements the negative profile loglikelihood function for the the parameter of interest **lambda** considering the conformable derivative.

**Input:**
- `VecPar124`: A vector containing the parameters **x0**, **K**, and **alpha**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.
- `VecPar3`: A numeric value for the parameter of interest, **lambda**.

 # 15. FunctionNegativeLogProfileLikFDElambda_Caputo.R
**Description:** Implements the negative profile loglikelihood function for the the parameter of interest **lambda** considering the Caputo derivative.

**Input:**
- `VecPar124`: A vector containing the parameters **x0**, **K**, and **alpha**.
- `VecTimeObs`: A vector of time points associated with each observation.
- `VecObs`: A vector of observed data.
- `h`: The step size for the computation.
- `m0`: The total number of points to generate.
- `VecPar3`: A numeric value for the parameter of interest, **lambda**.

# 16. FunctionProfileLCI.R
**Description:**  
Computes the profile likelihood confidence intervals for a specified level.

**Inputs:**  
- `Relative`: A vector containing the relative likelihood values.  
- `ValoresParametro`: A vector containing the parameter values associated with the relative likelihood.  
- `Nivel`: The desired confidence level for the likelihood interval, denoted as **c**.

  
