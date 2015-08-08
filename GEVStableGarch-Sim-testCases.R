
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
# TEST CASE FOF FUNCTION:               SPECIFICATION:
#  gsSim                        Useful for finding erros. Please add the expected 
#  						                  output of the test case whenever possible. 
#							                  This will help to trace erros during debugging.
################################################################################
# required packages to run test cases.
library("fGarch")
library("fExtremes")
library("stabledist")
library("skewt")
library("Rsolnp")



# ------------------------------------------------------------------------------
# Compare with the garchSim function using a presample matrix
# ------------------------------------------------------------------------------


# simulate GARCH(1,1) with conditional stable norm and compare with garchSim from fGarch.
# Note that the garchSim function from package fGarch has a somewhat strange instruction 
# < data = data[-(1:n.start), ] > works unexpectedly when n.start = 0. In these cases 
# we obtain a time series with size (n-1). Therefore, in order to compare our modified
# function gsSpec we need to put n.start > 0. 
presampleMatrix = matrix(c(0.1,0.5,0),1,3)
spec1 <- gsSpec(model = list(omega = 0.1, alpha = c(0.05), beta = c(0.02)), rseed = 1001, 
                     presample = presampleMatrix,cond.dist = "std")  
sim1 <- gsSim(spec = spec1, n = 5, n.start = 1)

spec2 <- garchSpec(model = list(omega = 0.1, alpha = c(0.05), beta = c(0.02)), 
                   presample = presampleMatrix,cond.dist = "std",rseed = 1001)  
sim2 <- .garchSim2(spec2, n = 5, n.start = 1, extended = TRUE)
# These means must be equal
mean(sim1)
mean(sim2)


# Test bigger matrix as a presample
presampleMatrix = matrix(c(0.1,-0.4,0.3,0.5,0,4),2,3)
spec1 <- gsSpec(model = list(omega = 0.1, alpha = c(0.05), beta = c(0.02)), rseed = 1001, 
                presample = presampleMatrix,cond.dist = "std")  
sim1 <- gsSim(spec = spec1, n = 5, n.start = 1)

spec2 <- garchSpec(model = list(omega = 0.1, alpha = c(0.05), beta = c(0.02)), 
                   presample = presampleMatrix,cond.dist = "std",rseed = 1001)  
sim2 <- garchSim2(spec2, n = 5, n.start = 1, extended = TRUE)
# These means must be equal
mean(sim1)
mean(sim2)


presampleMatrix = matrix(c(0.1,-0.4,-0.4,0.1,0.3,0.5,0,4,-4),3,3)
spec1 <- gsSpec(model = list(ar = c(-0.3,0.3), mu = -0.3, omega = 0.1, gamma = c(-0.3), alpha = c(0.05), 
               beta = c(0.02,0.3,0.3), delta = 0.5), rseed = 1001, 
               presample = presampleMatrix,cond.dist = "std")  
sim1 <- gsSim(spec = spec1, n = 5, n.start = 1)

spec2 <- garchSpec(model = list(ar = c(-0.3,0.3), mu = -0.3, omega = 0.1, gamma = c(-0.3), alpha = c(0.05), 
                                beta = c(0.02,0.3,0.3), delta = 0.5), 
                   presample = presampleMatrix,cond.dist = "std",rseed = 1001)  
sim2 <- garchSim2(spec2, n = 5, n.start = 1, extended = TRUE)
# These means must be equal
mean(sim1)
mean(sim2)

# ------------------------------------------------------------------------------
# Simulate GARCH(1,1) with several conditional distributions using 
# a set of parameters that makes the innovation and the simulated series
# be the same.
# ------------------------------------------------------------------------------
presample.matrix = matrix(c(0,1,0),1,3)
cond.dist.list = c("norm", "std", "sstd", "ged")
for( i in 1:length(cond.dist.list) )
{
    spec1 <- gsSpec(model = list(omega = 1, alpha = 0, beta = 0), rseed = 1001, 
                    presample = presample.matrix,cond.dist = cond.dist.list[i])  
    sim1 <- gsSim(spec = spec1, n = 5, n.start = 1)
    
    spec2 <- garchSpec(model = list(omega = 1, alpha = 0, beta = 0), rseed = 1001, 
                       presample = presample.matrix,cond.dist = cond.dist.list[i])
    sim2 <- garchSim2(spec2, n = 5, n.start = 1, extended = TRUE)
    
    if (sum( sim1 == sim2 ) != 15)
        stop("sim1 != sim2") 
}


# ------------------------------------------------------------------------------
# Testing the formula object 
------------------------------------------------------------------------------

# Expect a aparch model
spec = gsSpec(model = list(omega = 1, alpha = c(0.2), beta = c(0.5), 
	 skew = -0.5, shape = 1.5, delta = 1), cond.dist = "norm")
spec = gsSpec(model = list(omega = 1, alpha = c(0.2), beta = c(0.5), 
	 skew = -0.5, shape = 1.5, gamma = c(0.1)), cond.dist = "stableS1")
spec = gsSpec(model = list(omega = 1, alpha = c(0.2), beta = c(0.5), 
	 skew = -0.5, shape = 1.5, gamma = c(0.1)), cond.dist = "stableS1")
spec = gsSpec(model = list(omega = 1, alpha = c(0.2), beta = c(0.5), 
	 skew = -0.5, shape = 1.5, delta = 2), cond.dist = "stableS0")

# Expect a garch model
spec = gsSpec(model = list(omega = 1, alpha = c(0.2), beta = c(0.5), 
	 skew = -0.5, shape = 1.5, gamma = c(-0)), cond.dist = "stableS1")
spec = gsSpec(model = list(omega = 1, alpha = c(0.2), beta = c(0.5), 
	 skew = -0.5, shape = 1.5), cond.dist = "stableS1")
               