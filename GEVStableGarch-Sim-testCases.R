
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
#  GSgarch.Sim               Useful for finding erros. Please add the expected 
#  						 output of the test case whenever possible. 
#							 This will help to trace erros during debugging.
################################################################################
# required packages to run test cases.
library("fGarch")
library("fExtremes")
library("stabledist")
library("skewt")
library("Rsolnp")

# Simulate AR(1)-GARCH(1,1) with conditional GEV distr.
spec <- GSgarchSpec(model = list(ar = 0.1, alpha = 0.3, beta = 0.2, delta = 2), 
                    presample = NULL,cond.dist = "gev",rseed = NULL)  
sim <- GSgarch.Sim(spec, n = 1000, n.start = 0)
plot(sim)

# Simulate GARCH(1,1) with conditional stable distr.
spec <- GSgarchSpec(model = list(alpha = 0.05, beta = 0.01, omega = 0.01, delta = 2,shape = 1.2), 
                    presample = NULL,cond.dist = "stable",rseed = NULL)  
sim <- GSgarch.Sim(spec, n = 100, n.start = 0)
plot(sim)

# Simulate with small sample size: Expect error about size of n.
# spec <- GSgarchSpec(model = list(alpha = 0.05, beta = 0.01, omega = 0.01, delta = 2,shape = 1.2), 
#presample = NULL,cond.dist = "stable",rseed = NULL)  
#sim <- GSgarch.Sim(spec, n = 1, n.start = 0)

# simulate pure ARMA(1,1) with conditional stable distribution
spec <- GSgarchSpec(model = list(ar = 0.05, ma = 0.01,shape = 1.5,skew = -0.5), 
                    presample = NULL,cond.dist = "stable",rseed = NULL)  
sim <- GSgarch.Sim(spec, n = 100, n.start = 0)
plot(sim)

# simulate pure ARMA(2,2) with conditional stable distribution
spec <- GSgarchSpec(model = list(ar = c(0.05,0.3), ma = c(0.01,0.02),shape = 0.5,skew = -0.5), 
                    presample = NULL,cond.dist = "nor",rseed = NULL)  
sim <- GSgarch.Sim(spec, n = 100, n.start = 0)
plot(sim)

# simulate GARCH(1,1) with conditional stable norm and compare with garchSim from fGarch.
# Note that the garchSim function from package fGarch has a somewhat strange instruction 
# < data = data[-(1:n.start), ] > works unexpectedly when n.start = 0. In these cases 
# we obtain a time series with size (n-1). Therefore, in order to compare our modified
# function GSgarchSpec we need to put n.start > 0. 
presampleMatrix = matrix(c(0.1,0.5,0),1,3)
spec1 <- GSgarchSpec(model = list(omega = 0.1, alpha = c(0.05), beta = c(0.02)), 
                     presample = presampleMatrix,cond.dist = "std",rseed = 1001)  
sim1 <- GSgarch.Sim(spec1, n = 100, n.start = 1)

spec2 <- garchSpec(model = list(omega = 0.1, alpha = c(0.05), beta = c(0.02)), 
                   presample = presampleMatrix,cond.dist = "std",rseed = 1001)  
sim2 <- garchSim(spec2, n = 100, n.start = 1, extended = TRUE)
# These means must be equal
mean(sim1)
mean(sim2)                    


                    