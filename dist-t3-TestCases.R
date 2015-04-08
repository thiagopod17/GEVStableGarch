
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
#  pt3
#  dt3
################################################################################


# ------------------------------------------------------------------------------
# Test Cases for the density function dt3
# ------------------------------------------------------------------------------

# Test invalid input parameters 
dt3(rnorm(1000), mean = 0, sd = -1, nu = 2, d = 3, xi = 1)
dt3(rnorm(1000), mean = 0, sd = 3, nu = 0, d = -1, xi = 1)
dt3(rnorm(1000), mean = 0, sd = 3, nu = 1, d = -0.1, xi = 1)
dt3(rnorm(1000), mean = 0, sd = 3, nu = 0.001, d = 0.1, xi = 1)


# Test efficiency of implementation
n = 1e7; mean = rnorm(1,mean = 10); sd = runif(1,0,10); nu = runif(1,0,10); 
d = runif(1,0,10); xi = runif(1,0,10); x = rnorm(n)
system.time(
  dt3(x,  mean = mean, sd = sd, nu = nu, d = d, xi = xi)
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
  integrationValues[i] = as.numeric(integrate(dt3 , lower = -Inf, upper = Inf, 
                         mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                         d = dValues[i], xi = xiValues[i])[1])
}
error = 100*abs((integrationValues - trueValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(meanValues,sdValues,nuValues,dValues,xiValues,integrationValues)[which(error > 0.1),]
summary(error)



# ------------------------------------------------------------------------------
# Test Cases for the Distribution function pt3
# ------------------------------------------------------------------------------


# Plot some graphs for the distribution function:
n = 1e4; mean = rnorm(1,mean = 10); sd = runif(1,0,10); nu = runif(1,0,10); 
d = runif(1,0.1,10); xi = runif(1,0,1); x = sort(runif(n,-50,50))
y = pt3(x = x, mean = mean, sd = sd, nu = nu, d = d, xi = xi)
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
pt3LargeX = rep(NA,n)
pt3SmallX = rep(NA,n)
for(i in 1:n)
{
    pt3SmallX[i] = pt3 (smallX, mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                      d = dValues[i], xi = xiValues[i])
    pt3LargeX[i] = pt3 (largeX, mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
    d = dValues[i], xi = xiValues[i])  
}
errorSmall = 100*abs((pt3SmallX - 0))
errorLarge = 100*abs((pt3LargeX - 1))
summary(errorSmall)
summary(errorLarge)
cbind(meanValues,sdValues,nuValues,dValues,xiValues,pt3SmallX,pt3LargeX)[which(errorSmall > 0.05),]

# Test the computation of probability using the distribution function: 
# We will test the following property: integral( pt3, a, b) = pt3(b) - pt3(a)
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
  integrationValues[i] = as.numeric(integrate(dt3 , lower = aValues[i], upper = bValues[i], 
                                              mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                                              d = dValues[i], xi = xiValues[i])[1])
  trueValues[i] = pt3 (bValues[i], mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                      d = dValues[i], xi = xiValues[i]) - 
                  pt3 (aValues[i], mean = meanValues[i], sd = sdValues[i], nu = nuValues[i],
                  d = dValues[i], xi = xiValues[i])
}
error = 100*abs((integrationValues - trueValues))
summary(error)
plot(error, main = "Absolute error (%)", type = "l", col = "red")
cbind(meanValues,sdValues,nuValues,dValues,xiValues,aValues,bValues,trueValues,integrationValues)[which(error > 0.1),]