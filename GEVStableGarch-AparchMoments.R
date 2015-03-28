
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



####################################################################################
#  FUNCTION:                        DESCRIPTION:
#
#  Stationarity.Condition.Aparch    Evaluate the moments of  the type
#                                   E(z - gm|z|)^delta for several conditional
#                                   distributions.
####################################################################################

Stationarity.Condition.Aparch <-
  function (model = list(), 
            formula,
            cond.dist = c("gev","stable","norm", "std", "sstd"))
  {
    
    # Description: 
    #   A function that evaluate stationarity conditions 
    #   to guide parameter estimation in GARCH/APARCH models
    
    # Arguments:
    #   formula - an object returned by function .getFormula
    #   parm - a list with the model parameters as entries
    #     alpha - a vector of autoregressive coefficients
    #       of length p for the GARCH/APARCH specification,
    #     gm - a vector of leverage coefficients of
    #       length p for the APARCH specification,
    #     beta - a vector of moving average coefficients of
    #       length q for the GARCH/APARCH specification,
    #     delta - the exponent value used in the variance equation.
    #     skew - a numeric value listing the distributional
    #        skewness parameter.
    #     shape - a numeric value listing the distributional
    #        shape parameter.
    #   cond.dist - a character string naming the distribution
    #       function.           
    # make sure we have a garch or aparch with p > 0.
    
    
    # FUNCTION:
    
    # get parameters
    alpha <- model$alpha
    beta <- model$beta
    delta <- model$delta
    gm <- model$gm
    skew <- model$skew
    shape <- model$shape
    
    # Conditional distribution
    cond.dist = match.arg(cond.dist)
    
    if(length(alpha) == 0  || length(alpha) != length(gm) || length(delta) != 1)
      stop("Failed to verify conditions:
           if(length(alpha) == 0  || length(alpha) != length(gm) || length(delta) != 1)")
    
    # We must have a model with a non zero garch/aparch order
    if(formula$formula.order[3] == 0)
      stop("Invalid model")
    
    # garch model
    if(formula$isAPARCH == FALSE) 
      return(sum(alpha) + sum(beta))
    
    # aparch model
    kappa = rep(0,length(alpha))
    if(cond.dist == "norm")
      kappa <- norm.moment.aparch(delta = delta, gm = gm)
    if(cond.dist == "stable")
    {
      kappa <- stable.moment.aparch(shape = shape, skew = skew, 
                                    delta = delta, gm = gm)
    }
    # Distribution not implemented yet
    if( any(c("gev","std","sstd") == cond.dist) )   
      stop ("Stationarity.Condition.Aparch can not handle 
            this conditional distribution yet.")
    
    return(sum(kappa*alpha) + sum(beta))    
  }



# ------------------------------------------------------------------------------



norm.moment.aparch <- function(delta = 1.2, gm = 0)
{
  # Description:
  #   Returns the following Expectation for a standard normal distribution
  #   E[ (|z|-gm*z)^delta.
  #   Reference: A long memory property of stock market returns and a 
  #   new model. Ding, Granger and Engle (1993), Appendix B. 
  
  # Error treatment of input parameters
  if( (abs(gm) >= 1) || (delta <= 0) )
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true
         (abs(gm) >= 1) || (delta <= 0)")
  1/sqrt(2*pi)*( (1+gm)^delta + (1-gm)^delta )*
    2^((delta-1)/2)*gamma((delta+1)/2)
}



# ------------------------------------------------------------------------------



stable.moment.aparch <- function(shape = 1.5, skew = 0.5, delta = 1.2, gm = 0)
{
  # Description: 
  #   Returns the following Expectation for a standard normal distribution
  #   E[ (|z|-gm*z)^delta.
  #   This formula uses S(alpha,skew,1,0;pm) where pm = 2.
  #   Reference: The GEVStableGarch papper on JSS. 
  
  # Error treatment of input parameters
  if( (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) || 
        (delta <= 1) || (delta >= shape))
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true.
         (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || (abs(gm) >= 1) || 
         (delta <= 1) || (delta >= shape)
         Note: This function do not accept the normal case, i.e, alpha = 2")
  
  # Calculate the expression
  sigma.til <- (1 + (skew*tan(shape*pi/2))^2)^(1/2/shape)
  k.shape <- shape - 2
  skew.til <- 2/pi/(shape-2)*atan(skew*tan((shape-2)*pi/2))
  g1 <- gamma((1/2 + skew.til*k.shape/2/shape)*(-delta))
  g2 <- gamma(1/2 - skew.til*k.shape/2/shape + (1/2 + skew.til*k.shape/2/shape))
  g3 <- gamma((1/2 - skew.til*k.shape/2/shape)*(-delta))
  g4 <- gamma(1/2 + skew.til*k.shape/2/shape + (1/2 - skew.til*k.shape/2/shape))
  kappa <- 1/shape/sigma.til*sigma.til^(delta+1)*gamma(delta+1)*gamma(-delta/shape)*
    1/g1/g2*((1+gm)^delta + 1/g3/g4*(1-gm)^delta)
  return(kappa)
}



# ------------------------------------------------------------------------------
stable.simmetric.moment.garch <- function(shape = 1.5)
{
  # See Mittinik et al. (1995) for the definition of 
  # the stable garch model with conditional simmetric stable distribution
  # This formula uses S(alpha,skew,1,0;pm) where pm = 0 or 1. 
  if( (shape <= 1) || (shape > 2) )
    stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
         The following conditions cannot be true.
         (shape <= 0) || (shape > 2)")
  gamma(1 - 1/shape)*2/pi
}

# ------------------------------------------------------------------------------
stable.moment.power.garch <- function(shape = 1.5, skew = 0.5, delta = 1.2, gm = 0)
{
    # See Mittnik et al. (2002) - Stationarity of stable power-GARCH processes
    # This formula uses S(alpha,skew,1,0;pm) where pm = 0 or 1. 
    # Error treatment of input parameters
    if( (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || 
        (abs(gm) >= 1) || (delta >= shape))
      stop("Invalid parameters to calculate the expression E[ (|z|-gm*z)^delta ].
           The following conditions cannot be true.
            (shape <= 1) || (shape > 2) || ( abs(skew) > 1) || 
            (abs(gm) >= 1) || (delta >= shape)
           Note: This function do not accept the normal case, i.e, alpha = 2")
    
    
    # k(delta): At this point we have delta > 1
    if(delta == 1)
        k = pi/2
    else
        k = gamma(1-delta)*cos(pi*delta/2)
    # tal(,skew)
    tau = skew*tan(shape*pi/2)
    # lambda(shape,skew,delta): eq (10) Mittnik et al. (2002)
    lambda = k^(-1)*gamma(1-delta/shape)*(1 + tau^2)^(delta/2/shape)*cos(delta/shape*atan(tau))
    lambda
}


# ------------------------------------------------------------------------------



rest <- function(parm)
{
  # Description: Old and incomplete implementation of restriction of stationarity
  #   for aparch models
  mu <- parm[1];
  a <- parm[(1+1):(2+m-1)]; b <- parm[(1+m+1):(2+m+n-1)]
  omega <- parm[1+m+n+1]; alpha <- parm[(2+m+n+1):(3+m+n+p-1)]
  gm <- parm[(2+m+n+p+1):(3+m+n+p+p-1)]
  beta <- parm[(2+m+n+2*p+1):(3+m+n+2*p+q-1)]
  delta <- parm[2+m+n+2*p+q+1]; 
  skew <- parm[3+m+n+2*p+q+1]; shape <- parm[4+m+n+2*p+q+1];
  if(cond.dist == "stable")
  {
    tau <- skew*tan(shape*pi/2)
    kdelta <- pi/2
    if(abs(delta-1) > GStol) kdelta <- gamma(1 - delta)*cos(pi*delta/2)
    lamb <- kdelta^(-1)*gamma(1 - delta/shape)*(1 + tau^2)^(delta/2/shape)*
      cos(delta/shape*atan(tau))
    return(lamb*sum(alpha) + sum(beta))
  }
  return(sum(alpha) + sum(beta))
}



################################################################################