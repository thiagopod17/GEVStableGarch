
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
#  pGAt
#  dGAt
#  qGAt
#  rGAt
################################################################################


# ------------------------------------------------------------------------------
# Test Cases for the density function dGAt
# ------------------------------------------------------------------------------

# Test invalid input parameters 
dGAt(rnorm(1000), mean = 0, sd = -1, nu = 2, d = 3, xi = 1)
dGAt(rnorm(1000), mean = 0, sd = 3, nu = 0, d = -1, xi = 1)
dGAt(rnorm(1000), mean = 0, sd = 3, nu = 1, d = -0.1, xi = 1)
dGAt(rnorm(1000), mean = 0, sd = 3, nu = 0.001, d = 0.1, xi = 1)


# Test efficiency of implementation and compare with the dnorm function.
n = 1e7; mean = rnorm(1,mean = 10); sd = runif(1,0,10); nu = runif(1,0,10); 
d = runif(1,0,10); xi = runif(1,0,10); x = rnorm(n)
system.time(
  dGAt(x,  mean = mean, sd = sd, nu = nu, d = d, xi = xi)
  )
system.time(
  rnorm(n)
)

# Integrate the density and check if we get int(f,-infty,inifty) = 1
# 20150405 - Errors around 0.01% with maximum error equal to 0.12% of the true value.
# Notes: The integrate function fails to get the right value when d/sd is greater, 
# usually > 8.
# Indeed, small values of 'd' causes the integrate function to report errors or 
# wrong values. 
n = 10000
Smax = 10
Smin = 0.05
#meanValues = runif(n,-Smax,Smax)
#sdValues = runif(n,1,Smax)
meanValues = rep(0,n)
sdValues = rep(1,n)
nuValues = runif(n,Smin,Smax)
dValues = runif(n,0.3,Smax)
xiValues = runif(n,Smin,Smax)
integrationValues = rep(NA,n)
trueValues = 1
for(i in 1:n)
{
  integrationValues[i] = as.numeric(integrate(dGAt , lower = -Inf, upper = Inf, 
                         mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                         d = dValues[i], xi = xiValues[i])[1])
}
error = 100*abs((integrationValues - trueValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(meanValues,sdValues,nuValues,dValues,xiValues,integrationValues)[which(error > 0.1),]
summary(error)



# ------------------------------------------------------------------------------
# Test Cases for the Distribution function pGAt
# ------------------------------------------------------------------------------


# Plot some graphs for the distribution function:
n = 1e4; mean = rnorm(1,mean = 10); sd = runif(1,0,10); nu = runif(1,0,10); 
d = runif(1,0.1,10); xi = runif(1,0,1); x = sort(runif(n,-50,50))
y = pGAt(x = x, mean = mean, sd = sd, nu = nu, d = d, xi = xi)
plot(x,y, type = "l")

# Boundary tests:
# Tests the behaviour of the function when x approaches -Inf and +Inf
# 20150406 - These tests are valuable since they can help you in finding erros. 
# Change the largeX and smallX values to investigate the function behaviour. Since the domain
# of the distribution can be very sparse we may get P(X < smallX) = 0.1 even when smallX = 1e10.
# In these cases we need to test the function with smaller values, such as -1e-100 or sometimes
# -1e-400 to get values next to zero.
largeX = 1e10
smallX = -1e10
n = 10^5
Smax = 10
Smin = 0.05
meanValues = runif(n,-Smax,Smax)
sdValues = runif(n,1,Smax)
nuValues = runif(n,Smin,Smax)
dValues = runif(n,0.3,Smax)
xiValues = runif(n,Smin,1)
pGAtLargeX = rep(NA,n)
pGAtSmallX = rep(NA,n)
for(i in 1:n)
{
    pGAtSmallX[i] = pGAt (smallX, mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                      d = dValues[i], xi = xiValues[i])
    pGAtLargeX[i] = pGAt (largeX, mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
    d = dValues[i], xi = xiValues[i])  
}
errorSmall = 100*abs((pGAtSmallX - 0))
errorLarge = 100*abs((pGAtLargeX - 1))
summary(errorSmall)
summary(errorLarge)
cbind(meanValues,sdValues,nuValues,dValues,xiValues,pGAtSmallX,pGAtLargeX)[which(errorSmall > 0.05),]

# Test the computation of probability using the distribution function: 
# We will test the following property: integral( pGAt, a, b) = pGAt(b) - pGAt(a)
# Note: small values of 'd' may lead to numerical problems when computing the 
# probability with the 'integrate' routine from R.
# Indeed, when the density is very concentrated on some interval the integration routine
# may return a small (incorrect) value for the probability.
# 20150406 - Absolute errors around 2% when simulating 10^5 points. 
n = 10^5
Smax = 10
Smin = 0.05
rangeProbabilityValues = 1e1
aValues = rnorm(n,sd = rangeProbabilityValues)
bValues = aValues + abs(rnorm(n,sd = rangeProbabilityValues))
meanValues = runif(n,-Smax,Smax)
sdValues = runif(n,1,Smax)
nuValues = runif(n,Smin,Smax)
dValues = runif(n,0.3,Smax)
xiValues = runif(n,Smin,Smax)
integrationValues = rep(NA,n)
trueValues = rep(NA,n)
for(i in 1:n)
{
  integrationValues[i] = as.numeric(integrate(dGAt , lower = aValues[i], upper = bValues[i], 
                                              mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                                              d = dValues[i], xi = xiValues[i])[1])
  trueValues[i] = pGAt (bValues[i], mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                      d = dValues[i], xi = xiValues[i]) - 
                  pGAt (aValues[i], mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                  d = dValues[i], xi = xiValues[i])
}
error = 100*abs((integrationValues - trueValues))
summary(error)
plot(error, main = "Absolute error (%)", type = "l", col = "red")
cbind(meanValues,sdValues,nuValues,dValues,xiValues,aValues,bValues,trueValues,integrationValues)[which(error > 0.1),]


# ------------------------------------------------------------------------------
# Test Cases for the Quantile function qGAt
# ------------------------------------------------------------------------------


# Test the computation of quantiles and compare with the corresponding probability values
n = 100
nQuantiles = 5000
Smax = 10
Smin = 0.05
p = runif(nQuantiles,0,1)
meanValues = rnorm(n,0,Smax)
sdValues = runif(n,Smin,Smax)
nuValues = runif(n,Smin,Smax)
dValues = runif(n,0.3,Smax)
xiValues = runif(n,Smin,Smax)
quantileValues = matrix(NA,nrow = n, ncol = nQuantiles)
CalculatedValues = matrix(NA,nrow = n, ncol = nQuantiles)
error = matrix(NA,nrow = n, ncol = nQuantiles)
for(i in 1:n)
{
  quantileValues[i,] = qGAt(p,mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                                              d = dValues[i], xi = xiValues[i])
  CalculatedValues[i,] = pGAt(quantileValues[i,],mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                           d = dValues[i], xi = xiValues[i])
  error[i,] = 100*abs(CalculatedValues[i,] - p)/p
}

summary(as.vector(error))


# ------------------------------------------------------------------------------
# Test Cases for the random number generation function rGAt
# ------------------------------------------------------------------------------

# Plot of density and histogram
histPlot <- function(x, ...) {
  X = as.vector(x)
  H = hist(x = X, col = "steelblue",border = "white", nclass = 25, freq = FALSE)
  box()
  grid()
  abline(h = 0, col = "grey")
  mean = mean(X)
  sd = sd(X)
  xlim = range(H$breaks)
  s = seq(xlim[1], xlim[2], length = 201)
  lines(s, dGAt(s, ...), lwd = 2, col = "brown")
  abline(v = mean, lwd = 2, col = "orange")
  Text = paste("Mean:", signif(mean, 3))
  mtext(Text, side = 4, adj = 0, col = "darkgrey", cex = 0.7)
  rug(X, ticksize = 0.01, quiet = TRUE)
  invisible(s)
}
mean = 0; sd = 1; nu = 1; d = 5; xi = 1
random.GAt = rGAt(10000,mean = mean,sd = sd,nu = nu,d = d,xi = xi)
histPlot(random.GAt,
         mean = mean, sd = sd, nu = nu, d = d, xi = xi)

x = seq(-3,3,0.01)
plot(x , dGAt(x, mean = 0, nu = 1, d = 5, xi = 1))


# ------------------------------------------------------------------------------
# Example Scripts for the R help
# The examples were adapted from the examples given in the documentation 
# of the 'ged' density implemented inside package fGarch available on CRAN
# ------------------------------------------------------------------------------

# Simulate Random Values and compare with
# the empirical density and probability functions

# Configure plot and generate random values
par(mfrow = c(2, 2))
set.seed(1000)
r = rGAt(n = 1000)
plot(r, type = "l", main = "GAt Random Values", col = "steelblue")

# Plot empirical density and compare with true density:
hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue")
box()
x = seq(min(r), max(r), length = 201)
lines(x, dGAt(x), lwd = 2)

# Plot density function and compare with true df:
plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue",
     ylab = "Probability")
lines(x, pGAt(x), lwd = 2)

# Compute quantiles:
# Here we compute the quantiles corresponding to the probability points from 
# -10 to 10 and expect to obtain the same input sequence
round(qGAt(pGAt(x = seq(-10, 10, by = 0.5))), digits = 6)





################################################################################




