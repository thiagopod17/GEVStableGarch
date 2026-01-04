# USING DEVTOOLS

devtools::load_all() load your package (when run inside its folder)
devtools::test() to run all tests 

# Change log 

- allows only S0 parametrization for stable (cleaner code), and replaced dependence from stabledist by libstable4u.
- It allows exactly only stable, gev and gat distributions. Normal just for testing. 

# Tasks by priority

## Now
- Make testthat for garch(1,1) stable
- Build minimum test set that should go into pipeline (see tests folder):
    1) simulate and estimate using my package: model(m>=0,n>=0,p>=1,q>=0) + distributions (stableS0, GEV, GAT): test model combinations of (m,n,p,q) for m,n,q in [0,2] and p in [1,2] 
    2) testing extreme cases: simulate arma(1,1)-garch(1,1) with normal innovations. Estimate the same model using normal innovations, stable (alpha should be close to 2) and GAT (nu should be close to infinty) 
    1) garch(1,1) with normal dem2... with expected values
    2) test models and distribution
    2) simulate arma(1,1)-garch(1,1), fit and compare results. Do this for both stable, GEV and gat
    3) simulate arma(1,1)-garch(1,1) with normal and estimate with fGarch and rugarch. Same for stable conditional with alpha close to 2, and GAT with d=2, theta=1 should be t-student with v degrees of freedom. Do the oposite direction, simulate with fGarch and rugarch and estimate with mine. 


## Latter 
- Remove getFormula
- Remove Spec
- Simulation with libstable4u 
- Printing and computing hessian mixed inside gsFit function. gsFit is super big function. 
- My pdf Filtering Process for estimation (PDF DOC) is missing et = zt * ht in the equation
- Error message for computing std using hessian, not informative. users need mathematical reasons to investigage better the output of the function. 
- enforce stationarity using sqp.restriction algorithm must be tested with others datasets.






# Functions in this package

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

# Naming convention: 
user functions: gsFit, gsSelect, gsMomentAparch
variables: cond.dist, arma.order, garch.llh
constans: TOLG, TOLSTABLE, ARMA.ORDER
internal functions: .getStart, .getFormula
.

# Future modifications on package:

- Prediction methods using the results of Brockwell for stable prediction. See paper from Parameter "Estimation of ARMA Models with GARCH/APARCH Errors An R and SPlus Software Implementation"
  




