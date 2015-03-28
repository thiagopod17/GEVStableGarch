
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

# E(z^2) = Var(z) = 1, where z ~ N(0,1)
norm.moment.aparch(delta = 2, gm = 0)

# E(z^0) = 1, where z ~ N(0,1)
norm.moment.aparch(delta = 0.0000001, gm = 0)

# Comparison with simulated data
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
# Test Cases for function stable.simmetric.moment.garch
# ------------------------------------------------------------------------------

# E(|z|^1) = E(|x|^1) where x ~ N(0,sd = sqrt(2)) see Mittnik et al. (2002)
stable.simmetric.moment.garch(shape = 2)
sqrt(2)*norm.moment.aparch(delta = 1, gm = 0)

# Comparison with simulated data
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
                    stable.simmetric.moment.garch(shapeValues[i]))/
                    stable.simmetric.moment.garch(shapeValues[i]))
}
plot(error, ylab = "error In percentage")
summary(error)


# ------------------------------------------------------------------------------
# Test Cases for functions stable.moment.power.garch (Mittnik et al. (2002))
# ------------------------------------------------------------------------------
# E(|z|^1) = E(|x|^1) where x ~ N(0,sd = sqrt(2)) see Mittnik et al. (2002)
stable.moment.power.garch(shape = 2, skew = 0, delta = 1, gm = 0)
stable.simmetric.moment.garch(shape = 2)
sqrt(2)*norm.moment.aparch(delta = 1, gm = 0)

# Comparison with symmetric stable garch model
n = 10000
shapeValues = runif(n, 1, 2)
error = rep(NA,n)
for(i in 1:n)
{
  error[i] = abs((stable.moment.power.garch(shape=shapeValues[i],skew=0,delta=1,gm=0) - 
                    stable.simmetric.moment.garch(shape=shapeValues[i]) )/
                    stable.simmetric.moment.garch(shape=shapeValues[i]))
}
plot(error, ylab = "error In percentage", ylim = c(0,1e-15))
summary(error)
# ------------------------------------------------------------------------------
# Test Cases for functions stable.moment.aparch (GEVStableGarch papper)
# ------------------------------------------------------------------------------

# E(z^1) = E(x^1) where x ~ N(0,2) see Mittnik et al. (2002)
# CORRECT THIS FUNCTION BY STUDYING THE DIFFERENT PARAMETRIZATIONS ON THE 
# NOLAN BOOK CHAPTER 1.
stable.moment.aparch(shape = 1.6,skew = 0,delta = 1.5,gm = 0.5)
stable.moment.power.garch(shape = 1.6, skew = 0, delta = 1.5, gm = 0.5)
















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


