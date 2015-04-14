
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# TEST CASES FOF FUNCTIONS: 
#                
#  Stationarity.Condition.Aparch
#  norm.moment.aparch
#  stable.moment.aparch
################################################################################



# ------------------------------------------------------------------------------
# Test Cases for functions norm.moment.aparch
# ------------------------------------------------------------------------------

# E(z^2) = Var(z) = 1, where z ~ N(0,1) - OK
norm.moment.aparch(delta = 2, gm = 0)

# E(z^0) = 1, where z ~ N(0,1) - OK
norm.moment.aparch(delta = 0.0000001, gm = 0)

# Comparison with simulated data - OK
n = 100
nSimulations = 1000000
deltaValues = runif(n, 0, 10)
gmValues = runif(n, -1, 1)
error = rep(NA,n)
for(i in 1:n)
{
  z = rnorm(nSimulations)
  error[i] = abs((mean((abs(z)-gmValues[i]*z)^deltaValues[i]) - 
                    norm.moment.aparch(delta = deltaValues[i], gm = gmValues[i]))/
                   norm.moment.aparch(delta = deltaValues[i], gm = gmValues[i]))
}
plot(error, ylab = "error In percentage")
summary(error)




# ------------------------------------------------------------------------------
# Test Cases for functions gev.moment.aparch
# ------------------------------------------------------------------------------

# E(z^0) = 1, where z ~ GEV(0,1) - OK
gev.moment.aparch(delta = 0.0000001, gm = 0, shape = 5)

# Comparison with simulated data - OK
n = 10
nSimulations = 1000000
qsiValues = runif(n,-0.5,3)
deltaValues = rep(NA,n)
for(i in 1:n)
  deltaValues[i] = runif(1,0,1/qsiValues[i])
gmValues = runif(n, -1, 1)
error = rep(NA,n)
for(i in 1:n)
{
  z = rnorm(nSimulations)
  error[i] = abs((mean((abs(z)-gmValues[i]*z)^deltaValues[i]) - 
                    norm.moment.aparch(delta = deltaValues[i], gm = gmValues[i]))/
                   norm.moment.aparch(delta = deltaValues[i], gm = gmValues[i]))
}
plot(error, ylab = "error In percentage")
summary(error)





# ------------------------------------------------------------------------------
# Test Cases for functions std.moment.aparch (standard t-Student)
# ------------------------------------------------------------------------------


# Tests with the numerical integration computation - OK
n = 1000
shapeValues = runif(n, 2+0.01, 5)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,0,shapeValues[i]-0.09)
gmValues = runif(n,-1,1)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(TrueAparchMomentsWurtz(fun = "dstd",gm = gmValues[i], 
                                                    delta = deltaValues[i], nu = shapeValues[i])[1])
  functionValues[i] = std.moment.aparch(shape = shapeValues[i], delta = deltaValues[i],
                                        gm = gmValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(deltaValues,shapeValues,gmValues,trueValues,functionValues)
summary(error)



# ------------------------------------------------------------------------------
# Test Cases for functions ged.moment.aparch (standard GED distribution)
# ------------------------------------------------------------------------------

# Tests with the numerical integration computation - OK
n = 100
shapeValues = runif(n, 0.1, 5)
deltaValues = runif(n,0,10)
gmValues = runif(n,-1,1)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(TrueAparchMomentsWurtz(fun = "dged",gm = gmValues[i], 
                                                    delta = deltaValues[i], nu = shapeValues[i])[1])
  functionValues[i] = ged.moment.aparch(shape = shapeValues[i], delta = deltaValues[i],
                                        gm = gmValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(deltaValues,shapeValues,gmValues,trueValues,functionValues)
summary(error)


# ------------------------------------------------------------------------------
# Test Cases for functions sstd.moment.aparch (skew t-Student)
# ------------------------------------------------------------------------------
# We get numerical problems when delta/skew approaches 1.


# Tests with simulated data
# works well with errors around 0.001% of the true value.
# Sometimes the function TrueAparchMomentsWurtz fail because of the integration
# routine.


n = 1000
shapeValues = runif(n, 2, 15)
gmValues = runif(n,-1,1)
deltaValues = rep(NA,n)
skewValues = rep(NA,n)
for(i in 1:n)
{
  deltaValues[i] = runif(n=1,0.05,shapeValues[i]-0.2) 
  if(i %% 2 == 0)
    skewValues[i] = runif(n=1,0.015,deltaValues[i]-0.02)  
  else
    skewValues[i] = runif(n=1,deltaValues[i]+0.02,max(deltaValues[i]+20,20))
}
trueValues = rep(NA,n)
functionValues = rep(NA,n)
error = rep(NA,n)
for(i in 1:n)
{
  M1 = sqrt((shapeValues[i]-2)/pi)*gamma(shapeValues[i]/2)^(-1)*
    gamma((shapeValues[i]-1)/2)
  M2 = 1
  trueValues[i] = as.numeric(TrueAparchMomentsWurtz(fun = "dsstd", gm = gmValues[i], delta = deltaValues[i], 
                                                    nu = shapeValues[i], xi = skewValues[i], mean = (skewValues[i]-1/skewValues[i])*M1,
                                                    sd = sqrt((M2-M1^2)*(skewValues[i]^2+1/skewValues[i]^2)+2*M1^2-M2))[1])
  functionValues[i] = sstd.moment.aparch(shape = shapeValues[i], skew = skewValues[i], 
                                         delta = deltaValues[i], gm = gmValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(gmValues,skewValues,deltaValues,shapeValues,trueValues,functionValues)
summary(error)

# Compare the calculation of the moment with the G@RCH software
# For the dem2gbp-APARCH(1,1) the G@RCH reports: persistency = 0.9901654 but the true
# persistency is 0.9846524 as calculated with numerical integration on 
# the dsstd distribution

sstd.moment.aparch(shape = 4.211667, skew = exp(-0.100622), 
                   delta = 1.183602, gm = 0.149934)*0.135215 + 0.885756
.truePersistence(fun = "dsstd", alpha = 0.135215, beta = 0.885756,
                      gamma = 0.149934, delta = 1.183602, 
                       nu = 4.211667, xi = exp(-0.100622))

# For the dem2gbp-ARMA(1,1)-APARCH the G@RCH reports: persistency = 0.989947
# but the true persistency is 0.9845671.
sstd.moment.aparch(shape = 4.221416, skew = exp(-0.095899), 
                   delta = 1.202501, gm = 0.151121)*0.136845 + 0.884074
.truePersistence(fun = "dsstd", alpha = 0.136845, beta = 0.884074,
                 gamma = 0.151121, delta = 1.202501, 
                 nu = 4.221416, xi = exp(-0.095899))

# Conclusion: The G@RCH software is using the expression to calculate the 
# moments of skewed t-Student definied by Lambert and Laurent (2000, 2001),
# but the density defined in the garch model was reparameterized to be a 
# zero mean and unit variance.
# Another interesting conclusion is that if we work with the density 
# defined by Lambert and Laurent (without reparameterizing it to have 
# a zero mean and unit variance) we can also get the same garch estimated 
# parameters.



# ------------------------------------------------------------------------------
# Test Cases for functions t3.moment.aparch (standard t3 distribution)
# ------------------------------------------------------------------------------

# Tests with the numerical integration computation
# 20150413 - Test OK with error around 0.05% of the true value. 
# Sometimes the integrate function fails because of the integration interval.
# Note that the t3 density becames very picky when d approaches zero and thus,
# the integrate function will fail for these values. 
n = 2000
shapeValues = cbind(runif(n, 0.1, 5),runif(n, 0.5, 5))
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,0,shapeValues[i,1]*shapeValues[i,2]-0.09)
gmValues = runif(n,-1,1)
skewValues = runif(n,0,3)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(TrueAparchMomentsWurtz(fun = "dt3",gm = gmValues[i], 
                  delta = deltaValues[i], nu = shapeValues[i,1], d = shapeValues[i,2], 
                  xi = skewValues[i], lower = -Inf, upper = Inf)[1])
  functionValues[i] = t3.moment.aparch(shape = shapeValues[i,], delta = deltaValues[i],
                                        gm = gmValues[i], skew = skewValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(deltaValues,shapeValues,gmValues,trueValues,functionValues)[which(error > 1),]
summary(error)




# ------------------------------------------------------------------------------
# Test Cases for function stable.symmetric.moment.garch
# ------------------------------------------------------------------------------

# E(|z|^1) = E(|x|^1) where x ~ N(0,sd = sqrt(2)) see Mittnik et al. (2002) - OK
stable.symmetric.moment.garch(shape = 2)
sqrt(2)*norm.moment.aparch(delta = 1, gm = 0)

# Comparison with simulated data - OK
library(stabledist)
n = 100
nSimulations = 10000000
shapeValues = runif(n, 1.5, 2) # the mean of simulated values is more well behaved.
error = rep(NA,n)

for(i in 1:n)
{
  z = stabledist::rstable(n=nSimulations,alpha=shapeValues[i],beta = 0,
                          gm=1,delta=0,pm=1)
  error[i] = abs((mean(abs(z)) - 
                    stable.symmetric.moment.garch(shapeValues[i]))/
                   stable.symmetric.moment.garch(shapeValues[i]))
}
plot(error, ylab = "error In percentage")
summary(error)


# Comparison with the numerical integration routine
library(stabledist)
n = 20
shapeValues = runif(n, 1, 2)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
error = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(TrueAparchMomentsWurtz("dstable",gm = 0, delta = 1,
                                                    alpha = shapeValues[i], beta = 0,pm = 1)[1])
  functionValues[i] = stable.symmetric.moment.garch(shapeValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(shapeValues,trueValues,functionValues)
summary(error)







# ------------------------------------------------------------------------------
# Test Cases for functions stable.moment.power.garch (Mittnik et al. (2002))
# ------------------------------------------------------------------------------

# E(|z|^1) = E(|x|^1) where x ~ N(0,sd = sqrt(2)) see Mittnik et al. (2002) - OK
stable.moment.power.garch(shape = 2, skew = 0, delta = 1)
stable.symmetric.moment.garch(shape = 2)
sqrt(2)*norm.moment.aparch(delta = 1, gm = 0)

# Comparison with symmetric stable garch model - OK
n = 10000
shapeValues = runif(n, 1, 2)
error = rep(NA,n)
for(i in 1:n)
{
  error[i] = abs((stable.moment.power.garch(shape=shapeValues[i],skew=0,delta=1) - 
                    stable.symmetric.moment.garch(shape=shapeValues[i]) )/
                   stable.symmetric.moment.garch(shape=shapeValues[i]))
}
plot(error, ylab = "error In percentage", ylim = c(0,1e-15))
summary(error)

# Comparison with simulated data - OK
# Works really well for all values with errors around 0.002%!!!
# The trueValues are the estimated moments using numerical integration 
# from function TrueAparchMomentsWurtz addapted from function .truePersistence
# on the garch-Stats.R file.
n = 100
shapeValues = runif(n, 1+0.1, 2)
skewValues = runif(n,-1,1)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,1,shapeValues[i]-0.09)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(TrueAparchMomentsWurtz("dstable",gm = 0, delta = deltaValues[i],
                                                    alpha = shapeValues[i], beta = skewValues[i],pm = 1)[1])
  functionValues[i] = stable.moment.power.garch(shape = shapeValues[i], skew = skewValues[i],
                                                delta = deltaValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(skewValues,deltaValues,shapeValues,trueValues,functionValues)
summary(error)



# ------------------------------------------------------------------------------
# Test Cases for functions stable.moment.aparch (GEVStableGarch papper)
# ------------------------------------------------------------------------------

# E(z^1) = E(x^1), where x ~ N(0,2) see Mittnik et al. (2002) - OK
z = stabledist::rstable(n=10^7,alpha=1.999,beta = 0,
                        gm=1,delta=0,pm=0)
lambdaSim = mean((abs(z)-gm*z)^delta)
lambdaSim
stable.moment.aparch(shape = 1.999999,skew = skew,delta = delta,gm = gm)

# Comparison with stable.moment.power.garch function - OK
# from Mittnik et al. (2002)
n = 10000
shapeValues = runif(n, 1+0.1, 2)
skewValues = runif(n,-1,1)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,1,shapeValues[i]-0.09)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = stable.moment.power.garch(shape = shapeValues[i], skew = skewValues[i],
                                            delta = deltaValues[i])
  functionValues[i] = stable.moment.aparch(shape = shapeValues[i], skew = skewValues[i],
                                           delta = deltaValues[i], gm = 0)
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% GEVStableGarch eq. VS Mittnik)", type = "l", col = "red")
cbind(deltaValues,shapeValues,skewValues,trueValues,functionValues)
summary(error)

# Comparison with simulated data - OK
# Works really well for all values with errors around 0.005%!!!
n = 100
shapeValues = runif(n, 1+0.1, 2)
skewValues = runif(n,-1,1)
gmValues = runif(n,-1,1)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,1,shapeValues[i]-0.05)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(TrueAparchMomentsWurtz("dstable",gm = gmValues[i], 
                                                    delta = deltaValues[i],alpha = shapeValues[i], beta = skewValues[i],
                                                    pm = 1)[1])
  functionValues[i] = stable.moment.aparch(shape = shapeValues[i], skew = skewValues[i],
                                           delta = deltaValues[i], gm = gmValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(deltaValues,shapeValues,gmValues,skewValues,trueValues,functionValues)
summary(error)



# ------------------------------------------------------------------------------
# Test Cases for function stable.symmetric.moment.aparch (Diongue papper)
# ------------------------------------------------------------------------------

# Comparison with stable.moment.garch function 
# from GEVStableGarch papper on JSS - OK
# Error around 1e-9%!!!
n = 100000
shapeValues = runif(n, 1+0.1, 2)
gmValues = runif(n,-1,1)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,1,shapeValues[i]-0.09)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = stable.moment.aparch(shape = shapeValues[i], skew = 0,
                                       delta = deltaValues[i], gm = gmValues[i])
  functionValues[i] = stable.symmetric.moment.aparch(shape = shapeValues[i],
                                                     delta = deltaValues[i], gm = gmValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% GEVStableGarch eq. VS Mittnik)", type = "l", col = "red")
cbind(deltaValues,shapeValues,gmValues,trueValues,functionValues)
summary(error)




# ------------------------------------------------------------------------------
# Test Cases for function Stationarity.Condition.Aparch
# ------------------------------------------------------------------------------



# MA(2)-APARCH(1)-norm
spec <- GSgarchSpec(model = list(mu = 3,ma = c(1,2),alpha = 0.3,beta = 0.3, delta = 1), 
                    presample = NULL,cond.dist = c("norm"),rseed = 3)
# ARCH(1)-norm
spec <- GSgarchSpec(model = list(alpha = c(0.03), delta = 2), 
                    presample = NULL,cond.dist = c("norm"),rseed = 3)
# ARMA(2,3)-APARCH(2,2)-norm
spec <- GSgarchSpec(model = list(ar = c(1,2),ma = c(3,3,3), alpha = c(3,3),
                                 gm = c(0,0.4),beta = c(3,3),delta = 2), 
                    presample = NULL,cond.dist = c("norm"),rseed = 3)
# ARMA(2,3)-APARCH(2,2)-stable
spec <- GSgarchSpec(model = list(ar = c(0.1,0.04),ma = c(3,3,3), alpha = c(0.1,0.1),
                                 gm = c(0.3,0),beta = c(0.1,0.1),delta = 1.4, shape = 1.5, skew = 0), 
                    presample = NULL,cond.dist = c("stable"),rseed = 3)
Stationarity.Condition.Aparch(model = list(alpha = spec@model$alpha, beta = spec@model$beta, gm = spec@model$gm, 
                                           delta = spec@model$delta, skew = spec@model$skew, shape = spec@model$shape), 
                              formula = .getFormula(spec@formula), cond.dist = spec@distribution)

