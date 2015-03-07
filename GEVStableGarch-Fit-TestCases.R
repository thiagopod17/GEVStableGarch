
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
#  GSgarch.Fit                
#  Stationarity.Condition.Aparch
#  norm.moment.aparch
#  stable.moment.aparch
################################################################################



# ------------------------------------------------------------------------------
# Test Cases for functions GSgarch.Fit
# ------------------------------------------------------------------------------

############
# Testing different type of datasets as input
x <- c("asdf",rnorm(1000))
x <- c(rnorm(100),NA,rnorm(200))
x <- c(rnorm(100),-Inf)
x <- c(NULL)
x <- rnorm(100)
GSgarch.Fit(data = x , formula = ~arma(1,1)+garch(1,1),
            cond.dist = "norm", include.mean = TRUE, 
            algorithm = "sqp")

############
# Comparison with package GEVStableGarch from CRAN
# with dem2gbp[, 1] dataset
# Instructions: Clear the Workspace, 
# Load the packages and run 'model1', 2 and so on.
# Then load functions from the new version of the package 
# and run 'fit1', 2 and so on. Finally, compare the estimated
# parameters. 
library(fGarch)
library(GEVStableGarch)
data(dem2gbp)
x = dem2gbp[,1]

# garch(1,1)-norm-intercept
model1 = GSgarch.Fit(data = x, 0,0,1,1, cond.dist = "norm", 
                     intercept = TRUE, algorithm = "nlminb")

fit1 <- GSgarch.Fit(data = x , formula = ~garch(1,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb")

# garch(1,1)-std-intercept
model2 = GSgarch.Fit(data = x, 0,0,1,1, cond.dist = "std", 
                     intercept = TRUE, algorithm = "nlminb")

fit2 <- GSgarch.Fit(data = x , formula = ~garch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb")

# garch(1,1)-gev-intercept
model3 = GSgarch.Fit(data = x, 0,0,1,1, cond.dist = "gev", 
                     intercept = TRUE, algorithm = "nlminb")

fit3 <- GSgarch.Fit(data = x , formula = ~garch(1,1),
                    cond.dist = "gev", include.mean = TRUE, 
                    algorithm = "nlminb")

# arma(1,1)-aparch(1,1)-std-intercept
model4 = GSgarch.Fit(data = x, 1,1,1,1, cond.dist = "std", 
                     intercept = TRUE, algorithm = "nlminb",APARCH = TRUE)

fit4 <- GSgarch.Fit(data = x , formula = ~arma(1,1)+aparch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb")

model1@fit$par
fit1@fit$par
model2@fit$par
fit2@fit$par
model3@fit$par
fit3@fit$par
model4@fit$par
fit4@fit$par


############
# Comparison between GEVStableGarch and fGarch fit functions
# with dem2gbp[, 1] dataset

library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]

# garch(1,1)-norm-intercept
fit1 <- garchFit(data = x, formula = ~garch(1,1),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb")
model1 <- GSgarch.Fit(data = x , formula = ~garch(1,1),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "sqp")
fit1@fit$par-model1@fit$par
fit1@fit$llh
model1@fit$llh
# garch(1,1)-std-intercept
fit1 <- garchFit(data = x, formula = ~garch(1,1),
                      cond.dist = "std", include.mean = TRUE,
                      algorithm = "nlminb")
model1 <- GSgarch.Fit(data = x , formula = ~garch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb")
fit1@fit$par-model1@fit$par

# garch(2,2)-norm-intercept
fit1 <- garchFit(data = x, formula = ~garch(2,2),
                 cond.dist = "std", include.mean = TRUE,
                 algorithm = "nlminb")


model1 <- GSgarch.Fit(data = x , formula = ~garch(2,2),
                      cond.dist = "std", include.mean = TRUE, 
                      algorithm = "nlminb")
fit1@fit$par-model1@fit$par

# garch(1,0)-norm-intercept
fit1 <- garchFit(data = x, formula = ~garch(1,0),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb")


model1 <- GSgarch.Fit(data = x , formula = ~garch(1,0),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "nlminb")
fit1@fit$par-model1@fit$par

# aparch(1,1)-norm-intercept
fit1 <- garchFit(data = x, formula = ~aparch(1,1),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb")


model1 <- GSgarch.Fit(data = x , formula = ~aparch(1,1),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "nlminb", DEBUG = TRUE, control = list(trace = 3))
fit1@fit$par-model1@fit$par

# aparch(1,0)-norm-intercept
fit1 <- garchFit(data = x, formula = ~aparch(1,0),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb")


model1 <- GSgarch.Fit(data = x , formula = ~aparch(1,0),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "nlminb", DEBUG = FALSE)
fit1@fit$par-model1@fit$par

############
# Fitting ARMA-GARCH or ARMA-APARCH process
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]+100

# arma(1,1)-garch(1,1)-std-intercept
fit1 <- GSgarch.Fit(data = x, formula = ~arma(1,1)+garch(1,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- garchFit(data = x, formula = ~arma(1,1)+garch(1,0),
                    cond.dist = "norm", include.mean = TRUE)
model1@fit$par
fit1@fit$par
model1@fit$par-fit1@fit$par

# arma(5,0)-garch(1,0)-std-intercept
fit1 <- GSgarch.Fit(data = x, formula = ~arma(5,0)+garch(1,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- garchFit(data = x, formula = ~arma(5,0)+garch(1,0),
                   cond.dist = "norm", include.mean = TRUE)
model1@fit$par
fit1@fit$par
model1@fit$par-fit1@fit$par

# arma(0,5)-garch(1,0)-std-intercept
fit1 <- GSgarch.Fit(data = x, formula = ~arma(0,5)+garch(1,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- garchFit(data = x, formula = ~arma(0,5)+garch(1,0),
                   cond.dist = "norm", include.mean = TRUE)
model1@fit$par
fit1@fit$par
model1@fit$par-fit1@fit$par

# Fitting ARMA-GARCH or ARMA-APARCH process
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]

# arma(1,1)-garch(1,0)-norm-intercept
fit1 <- GSgarch.Fit(data = x, formula = ~arma(1,1)+garch(1,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- garchFit(data = x, formula = ~arma(1,1)+garch(1,0),
                   cond.dist = "norm", include.mean = TRUE)
model1@fit$par-fit1@fit$par


# arma(1,1)-aparch(1,1)-std-intercept
fit1 <- GSgarch.Fit(data = x, formula = ~arma(1,1)+aparch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "sqp")

model1 <- garchFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                   cond.dist = "std", include.mean = TRUE)
(model1@fit$par-fit1@fit$par)/fit1@fit$par



# arma(0,1)-garch(1,1)-norm-intercept
data(sp500dge)
x = 100*sp500dge[, 1]
fit1 <- GSgarch.Fit(data = x, formula = ~arma(0,1)+aparch(1,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "sqp")
model1 <- garchFit(data = x, formula = ~arma(0,1)+aparch(1,1),
                   cond.dist = "norm", include.mean = TRUE, 
                   algorithm = "nlminb")
model1@fit$par-fit1@fit$par

############
# Fitting models using different Algorithms ('nlminb' and 'sqp')

# arma(1,1)-aparch(1,1)-std
fit1 <- GSgarch.Fit(data = x, formula = ~arma(1,1)+aparch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- GSgarch.Fit(data = x, formula = ~arma(1,1)+aparch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "sqp")
abs((model1@fit$par-fit1@fit$par)/fit1@fit$par)

# arma(1,1)-std
fit1 <- GSgarch.Fit(data = x, formula = ~arma(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- GSgarch.Fit(data = x, formula = ~arma(1,1),
                      cond.dist = "std", include.mean = TRUE, 
                      algorithm = "sqp")

abs((model1model1@fit@par-fit1@fit$par)/fit1@fit$par)

# aparch(2,2)-std
fit1 <- GSgarch.Fit(data = x, formula = ~aparch(2,2),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- GSgarch.Fit(data = x, formula = ~aparch(2,2),
                      cond.dist = "std", include.mean = TRUE, 
                      algorithm = "sqp")

model2 <- garchFit(data = x, formula = ~aparch(2,2),
                      cond.dist = "std", include.mean = TRUE, 
                      algorithm = "nlminb")

abs((model1@fit$par-fit1@fit$par)/fit1@fit$par)
abs((model2@fit$par-model1@fit$par)/model1@fit$par)

# arma(0,1)-norm
fit1 <- GSgarch.Fit(data = x, formula = ~arma(0,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- GSgarch.Fit(data = x, formula = ~arma(0,1),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "sqp")

abs((model1@fit$par-fit1@fit$par)/fit1@fit$par)

# arma(1,0)-aparch(1,0)-norm
fit1 <- GSgarch.Fit(data = x, formula = ~arma(0,1)+aparch(1,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- GSgarch.Fit(data = x, formula = ~arma(0,1)+aparch(1,0),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "sqp")

abs((model1@fit$par-fit1@fit$par)/fit1@fit$par)

############
# Fitting pure ARMA process 
# Notes:
# The arma(m,0) or arma(0,n) are perfectly fitted by our algorithm. This 
# happens with both ARMA-GARCH models and ARMA only models.
# There is still a problem with the ARMA filter function for ARMA(1,1) and 
# other combinations.
# The "sigma" parameter is the square root of the sigma2 parameter estimated
# by function "arima".
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]+30
# arma(1,1)-norm-intercept-nlminb
m <- 4
n <- 0
fit1 <- GSgarch.Fit(data = x, formula = ~arma(4,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb", DEBUG = FALSE, control = list(trace = 3))
model1 <- arima(x, order = c(m, 0, n))
fit1@fit$par
model1$coef[c(m+n+1,1:(m+n))]
absoluteError <- abs((fit1@fit$par[1:(1+m+n)]-model1$coef[c(m+n+1,1:(m+n))])/fit1@fit$par[1:(1+m+n)])
absoluteError
fit1@fit$llh
model1$loglik
fit1 <- GSgarch.Fit(data = x, formula = ~arma(2,2),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb", DEBUG = FALSE)
model1 <- arima(x, order = c(2, 0, 2), include.mean = TRUE)
absoluteError <- abs((fit1$par-model1$coef[c(5,1:4)])/fit1$par)
absoluteError

############
# Testing the output object of class fGEVSTABLEGARCH

library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]
# garch(1,1)-norm-intercept
fit1 <- GSgarch.Fit(data = x, formula = ~garch(1,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb")
model1 <- garchFit(data = x, formula = ~garch(1,1),
                   cond.dist = "norm", include.mean = TRUE, 
                   algorithm = "nlminb")
model1@fit
fit1@fit

#########
# Studying the ARMA stationarity
# The arCheck function says TRUE if the ARMA model is stationary.
arCheck <- function(ar) {
  p <- max(which(c(1, -ar) != 0)) - 1
  if (!p) 
    return(TRUE)
  all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
}
arCheck(c(-0.61))
# antes do arCheck -0.02665051 -0.62383081  0.64553619 
# Depois do arCheck -0.02665051 -0.62383081  0.64553619 
# ------------------------------------------------------------------------------
# Test Cases for functions norm.moment.aparch and stable.moment.aparch
# ------------------------------------------------------------------------------



# E(z^2) = Var(z) = 1
norm.moment.aparch(delta = 2, gamma = 0)
stable.moment.aparch(2,0.5,3,0)



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
                                 gamma = c(0,0.4),beta = c(3,3),delta = 2), 
                    presample = NULL,cond.dist = c("norm"),rseed = 3)
# ARMA(2,3)-APARCH(2,2)-stable
spec <- GSgarchSpec(model = list(ar = c(0.1,0.04),ma = c(3,3,3), alpha = c(0.1,0.1),
                                 gamma = c(0.3,0),beta = c(0.1,0.1),delta = 1.4, shape = 1.5, skew = 0), 
                    presample = NULL,cond.dist = c("stable"),rseed = 3)
Stationarity.Condition.Aparch(model = list(alpha = spec@model$alpha, beta = spec@model$beta, gamma = spec@model$gamma, 
                                           delta = spec@model$delta, skew = spec@model$skew, shape = spec@model$shape), 
                              formula = .getFormula(spec@formula), cond.dist = spec@distribution)



################################################################################

