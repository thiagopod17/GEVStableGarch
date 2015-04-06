
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


n = 1e4; mean = rnorm(1,mean = 10); sd = runif(1,0,10); nu = runif(1,0,10); 
d = runif(1,0.1,10); xi = runif(1,0,1); x = sort(runif(n,-50,50))
y = pt3(x = x, mean = mean, sd = sd, nu = nu, d = d, xi = xi)
plot(x,y, type = "l")


