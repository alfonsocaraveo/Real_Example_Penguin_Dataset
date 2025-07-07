# Repository containing model selection procedures applied to a penguin dataset
Repository containing model selection procedures applied to the penguin dataset, 
based on the article by A. Chiaradia and I. C. Nisbet on plasticity in parental provisioning and 
chick growth in little penguins (Eudyptula minor) during years of high and low breeding success.

# Scripts Overview for the model selection analysis
This repository contains the following R scripts, each serving a specific purpose in the workflow:

# 1. .R
**Description:** Script for.

**Libraries Used:** 

- `lubridate`: For date and time manipulation.
- `mnormt`: For multivariate normal distribution functions.


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

