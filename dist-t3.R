
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
# FUNCTION:              t3 DISTRIBUTION PROPOSED IN PAOLELLA (1997)
#  pt3                   Probability function for the t3
#  dt3                   Density for the t3-distribution 
#  qt3                   Quantile function for the t3
#  rt3                   Random Number Generator for the t3
################################################################################


dt3 <- 
  function(x, mean = 0, sd = 1, nu = 2, d = 3, xi = 1, log = FALSE)
  {   
    # A function imlemented by Thiago Sousa
    
    # Description:
    #   Compute the density for the 
    #   so called t3-distribution defined in Paolella (1997).
    #   Reference: Paolella M 0886. Tail Estimation and Conditional Modeling 
    #   of Heteroscedastic Time!Series. PhD thesis.
    #   Institute of Statistics and Econometrics. Christian Albrechts University at Kiel
    #   Parameters: mean in R; sd > 0; nu > 0; d > 0; xi > 0; 
    
    # FUNCTION:
    
    # Params:
    if (length(mean) == 5) {
      xi = mean[5]
      d  = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }    
    
    # Error treatment of input parameters
    if(sd <= 0  || nu <= 0 || xi <= 0 || d <= 0)
      stop("Failed to verify condition:
           sd <= 0 || nu <= 0 || xi <= 0 || d <= 0")
    
    # Compute auxiliary variables:
    z = (x - mean ) / sd
    n = length(z)
    arg = z
    indexLessThanZero = which (z < 0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanZero)
    
    # Compute the density points according to their sign
    # all b coefficients are >= 0
    if(sizeIndex == 0) {
      arg = arg / xi  
    # all b coefficients are < 0
    } else if (sizeIndex == n) {
        arg = -arg * xi 
    # default case. we have both pos. and neg. values
    } else if (TRUE) { 
        arg[indexLessThanZero] = -arg[indexLessThanZero] * xi
        arg[-indexLessThanZero] = arg[-indexLessThanZero] / xi 
    }
    
    # Compute density points
    k = ( ( xi + 1/xi ) * 1/d * nu^(1/d) * beta (1/d,nu) )^(-1)
    result = ( k * (1 + (arg^d) / nu )^( -nu-1/d) ) / sd
    # Log:
    if(log) result = log(result)
    
    # Return Value
    result
  }





# ------------------------------------------------------------------------------

pt3 <- 
  function(x, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
  {   
    # A function imlemented by Thiago Sousa
    
    # Description:
    #   Compute the distribution for the 
    #   so called t3-distribution
    #   Parameters: mean in R; sd > 0; nu > 0; d > 0; xi > 0; 
    
    # FUNCTION:
    
    # Params:
    if (length(mean) == 5) {
      xi = mean[5]
      d  = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }    
    
    # Error treatment of input parameters
    if(sd <= 0  || nu <= 0 || xi <= 0 || d <= 0)
      stop("Failed to verify condition:
           sd <= 0 || nu <= 0 || xi <= 0 || d <= 0")
    
    # Define auxiliary functions
    L <- function (z, nu = nu, d = d, xi = xi)
    {
        nu / ( nu + (-z*xi)^d ) # z must be negative ( <= 0 ), but we do not check it here
    }
    U <- function (z, nu = nu, d = d, xi = xi)
    {
        pw = (z/xi)^d  # z must be negative ( <= 0 ), but we do not check it here
        if(pw == Inf)
            1
        else
            pw / ( nu + pw ) 
    } 

    # Compute auxiliary variables:
    z = (x - mean ) / sd
    n = length(z)
    arg = z
    indexLessThanZero = which (z <= 0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanZero)
    
    # Compute distribution points according to their sign
    if(sizeIndex == 0) {
        arg = 1/(1 + xi^2 ) + 1/(1 + xi^(-2) ) * 
            pbeta ( U (z = arg, nu = nu, d = d, xi = xi), 1/d, nu)  
    } else if (sizeIndex == n) {
        arg = 1/(1 + xi^2 ) * 
            pbeta ( L (z = arg, nu = nu, d = d, xi = xi), nu, 1/d)
    } else if (TRUE) { 
        arg[indexLessThanZero] = 1/(1 + xi^2 ) * 
            pbeta ( L (z = arg[indexLessThanZero], nu = nu, d = d, xi = xi), nu, 1/d)
        arg[-indexLessThanZero] = 1/(1 + xi^2 ) + 1/(1 + xi^(-2) ) * 
            pbeta ( U (z = arg[-indexLessThanZero], nu = nu, d = d, xi = xi), 1/d, nu) 
    }
    
    # Return Value
    arg
  }

# ------------------------------------------------------------------------------
# 
# 
# qged <- 
#   function(p, mean = 0, sd = 1, nu = 2)
#   {   
#     # A function implemented by Diethelm Wuertz
#     
#     # Description:
#     #   Compute the quantiles for the  
#     #   generalized error distribution.
#     
#     # FUNCTION:
#     
#     # Compute Quantiles:
#     lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
#     q = lambda * (2*qgamma((abs(2*p-1)), 1/nu))^(1/nu)
#     result = q*sign(2*p-1) * sd + mean
#     
#     # Return Value:
#     result
#   }
# 
# 
# # ------------------------------------------------------------------------------
# 
# 
# rged <-  
#   function(n, mean = 0, sd = 1, nu = 2)
#   {   
#     # A function implemented by Diethelm Wuertz
#     
#     # Description:
#     #   Generate GED random deviates. The function uses the 
#     #   method based on the transformation of a Gamma random 
#     #   variable.
#     
#     # FUNCTION:
#     
#     # Generate Random Deviates:
#     lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
#     # print(lambda)
#     r = rgamma(n, 1/nu)
#     z =  lambda * (2*r)^(1/nu) * sign(runif(n)-1/2)
#     result = z * sd + mean
#     
#     
#     # Return Value:
#     result
#   }


################################################################################

