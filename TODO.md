# TASKS BY PRIORITY

## build code 
- OK. CREATE MM WITH HIGHER ORDER STRUCTURE FROM TOPIC FUNCTIONS OF PACKAGE. THEN, CREATE A MAP FROM ALL FUNCTIONS AVAILABLE, TO SEE WHAT IS REALLY NECESSARY, THERE IS A LOT OF DUPLICATED CODE.

## testing 
- tests after building package: this R file should be merged, I dont want to have manual testing. Same holds for the testsCasesFiles, they seem manual.
- GEVStableGarch-armaDist.R: lots of replicated code, I wonder if is possible to make one code only setting a generic distribution. 

## build package
- advices or functions for building the package, maybe move to the other repository or to readme functions. They seem like instructions. 




# Functions of package

 .armaDist               Calculates the likelihood function for a vector of points (z)
    #   according to the specified distribution.

 .armaGarchDist          similar to .armaDist but density for z/hh

 .filterEQUATION.        filter data according to a specified model by EQUATION 

 Fit -  Fits ARMA-GARCH or ARMA-APARCH model. Returns an object of class GEVSTABLEGARCH 

 getFormula  - get model parameters from string, like, 'arma(1,1)-garch(2,2)'

 gsSelect (previous or related name = GSgarch.FitAIC) (returns model parameters with minimum AIC). It needs getOrder to get list of different parameters to fit and decide which one is the best. 

 getStart - starting parameters for estimation, relies on preliminary arima fitting 

otherUsefulCodes - for testing ? 

Sim - model simulation

Spec - specifies model and returns an instance of class GEVSTABLEGARCHSPEC

stationarityAparch - compute a value that should be < 1 for the model to be stationary.







# OLD TODO LIST, TAKE WHAT MAKES SENSE FROM HERE


#   Things to do on Package now
    - send package to CRAN
    - organize my files on repository
    - copy code from package to repository

# Improvements 
    - sqp.restriction algorithm must be tested with others datasets.

# Naming: 
      user functions: gsFit, gsSelect, gsMomentAparch
      variables: cond.dist, arma.order, garch.llh
      constans: TOLG, TOLSTABLE, ARMA.ORDER
      internal functions: .getStart, .getFormula


# Changes on function names:

gsGarchDist    .armaGarchDist
filter.Arma     .filterArma
filter.Aparch   .filterAparch
filter.Aparch.Forloop     .filterAparchForLoop
filter1.garch11Fit
gsGetOrder    .getOrder
gsGetStart     .getStart
changed variable name from 'gm' to 'gamma'











# Advices for Debugging: 

  - See the TOLG and TOLGSTABLE parameters in function gsGetStart. They were originally set to 
    1e-7 and 2e-2.





# Future modifications on package:

    - Implement the ARMA dist function for every distribution. 
    - Currently the code is working for ARMA(1,1) ARMA(p,1) ARMA(1,n) models with condtional normal. 
    - Find A More Efficiet Way To Calculate The GEV Aparch Moment Instead Of Using The Integration function.
    - Advices of professor Doctor Paolella.
    - Prediction methods using the results of Brockwell for stable prediction.
    - Advices professor Paolella. (email)
  




