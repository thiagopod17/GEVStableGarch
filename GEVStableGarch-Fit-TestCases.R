
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
# TEST CASES FOF FUNCTION:               SPECIFICATION:
#  GSgarch.Fit                Test on input parameters and on numerical 
#                             stability of parameter estimation 
#  						                							               
################################################################################

# Historical Notes.
# We will make the interface of this function more similar to the garchFit function
# from package fGarch to make its use easier. 
# The input parameters for this version (march/2015) are:

# 25 Feb, 2015.
# Parameter input modification
  # OK. algorithm: a string parameter that determines the algorithm used for maximum likelihood estimation.
  # OK. cond.dist: name of the conditional distribution, one of gev, stable, norm, std, sstd
  # OK. control: control parameters, the same as used for the functions from nlminb, and 'bfgs' and 'Nelder-Mead' from optim.
  # OK. data: The dataset to be estimated.
  # formula: formula object describing the mean and variance equation of the ARMA-GARCH/APARCH model.
  # OK. intercept: this flag determines if the parameter for the mean will be estimated or not
  # OK. print.Result (Padrao eh TRUE): A boolean variable specifying whether or not the user wants to print the results after the function calling.
  # OK. get.res: (NAO VAMOS TER MAIS ESSA VARIAVEL)
  # OK. GSstable.tol e GStol: (CONFIGURAR NO INICIO DA FUNCAO, NAO MAIS NECESSARIA AQUI)
  # APARCH: Vamos tirar pois usaremos da formula.
# 25 feb, 2015, right before commiting on Github
  # We saw that the estimated parameters from both GEVStableGarch package from CRAN
  # and from our current version are the same. On the other hand, we saw that the 
  # estimated parameters from macbook differ slightly from the windows version.
  # our goal now is to investigate the filtering process inside the GSGarch.Fit 
  # function to make the estimated parameters more similar to the ones from package
  # fGarch. 
# 27 feb, 2015, right before commiting to Github
  # Using garch11Fit function from Wurtz (2006) to estimate
  # pure garch(1,1) model with conditional normal distribution
  # This function estimate the
  # garch(1,1)-include.mean-norm-dem2gbp
  # The results are exactly the same as in the Code Snippet 2
  # presented in the papper Wurtz et al. (2006)
  # The function garch11Fit works better if start the conditional 
  # variance with 'var(x)'.
  # Mehoramos muito minha funcao quando para a estimacao do garch(1,1). Fiz
  # isso retirando o filtro do aparch e recolocando o filtro do wuertz, que funcionava
  # para o garch11.
  # When I put the filter from garch11Fit function inside my GSgarch.Fit the 
  # results were exactly the same. Therefore, thats is our starting pointing.
  # Now, I am almost done because my filter function for pure APARCH model is 
  # is matching exactly the filter function from garch11Fit. 
  # The next step is to test it considerably well and develop the other filtering 
  # for other process. This function needs to be documented a lot. Also, 
  # remember to take pictures of the matrix representation I did on paper to 
  # commit to github.


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

model1$par
fit1$par
model2$par
fit2$par
model3$par
fit3$par
model4$par
fit4$par


############
# Comparison between GEVStableGarch and fGarch fit functions
# with dem2gbp[, 1] dataset

library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]

# garch(1,1)-std-intercept
fit1 <- garchFit(data = x, formula = ~garch(1,1),
                      cond.dist = "std", include.mean = TRUE,
                      algorithm = "nlminb")
model1 <- GSgarch.Fit(data = x , formula = ~garch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb")
fit1@fit$matcoef[,1]-model1$matcoef[,1]

# garch(2,2)-norm-intercept
fit1 <- garchFit(data = x, formula = ~garch(2,2),
                 cond.dist = "std", include.mean = TRUE,
                 algorithm = "nlminb")


model1 <- GSgarch.Fit(data = x , formula = ~garch(2,2),
                      cond.dist = "std", include.mean = TRUE, 
                      algorithm = "nlminb")
fit1@fit$matcoef[,1]-model1$matcoef[,1]

# garch(1,0)-norm-intercept
fit1 <- garchFit(data = x, formula = ~garch(1,0),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb")


model1 <- GSgarch.Fit(data = x , formula = ~garch(1,0),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "nlminb")
fit1@fit$matcoef[,1]-model1$matcoef[,1]

# aparch(1,1)-norm-intercept
fit1 <- garchFit(data = x, formula = ~aparch(1,1),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb")


model1 <- GSgarch.Fit(data = x , formula = ~aparch(1,1),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "nlminb", DEBUG = TRUE, control = list(trace = 3))
fit1@fit$matcoef[,1]-model1$matcoef[,1]

# aparch(1,0)-norm-intercept
fit1 <- garchFit(data = x, formula = ~aparch(1,0),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb")


model1 <- GSgarch.Fit(data = x , formula = ~aparch(1,0),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "nlminb", DEBUG = FALSE)
fit1@fit$matcoef[,1]-model1$matcoef[,1]


############
# Fitting ARMA-GARCH or ARMA-APARCH process
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]

# arma(1,1)-garch(1,1)-std-intercept
fit1 <- GSgarch.Fit(data = x, formula = ~arma(1,1)+garch(1,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- garchFit(data = x, formula = ~arma(1,1)+garch(1,1),
                    cond.dist = "norm", include.mean = TRUE)
model1@fit$matcoef[,1]-fit1$matcoef[,1]


# arma(1,1)-aparch(1,1)-std-intercept
fit1 <- GSgarch.Fit(data = x, formula = ~arma(1,1)+aparch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb")

model1 <- garchFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                   cond.dist = "std", include.mean = TRUE)
model1@fit$matcoef[,1]-fit1$matcoef[,1]


# arma(0,1)-garch(1,1)-norm-intercept
data(sp500dge)
x = 100*sp500dge[, 1]
fit1 <- GSgarch.Fit(data = x, formula = ~arma(0,1)+aparch(1,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "sqp")
model1 <- garchFit(data = x, formula = ~arma(0,1)+aparch(1,1),
                   cond.dist = "norm", include.mean = TRUE, 
                   algorithm = "nlminb")
model1@fit$matcoef[,1]-fit1$matcoef[,1]

############
# Fitting pure ARMA process 
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]

# arma(1,1)-norm-intercept-nlminb
fit1 <- GSgarch.Fit(data = x, formula = ~arma(1,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb", DEBUG = FALSE)
model1 <- arima(x, order = c(1, 0, 1))
absoluteError <- abs((fit1$par-model1$coef[c(3,1,2)])/fit1$par)
absoluteError

# arma(2,2)-norm-intercept-nlminb
fit1 <- GSgarch.Fit(data = x, formula = ~arma(2,2),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb", DEBUG = FALSE)
model1 <- arima(x, order = c(2, 0, 2), include.mean = TRUE)
absoluteError <- abs((fit1$par-model1$coef[c(5,1:4)])/fit1$par)
absoluteError

################################################################################

