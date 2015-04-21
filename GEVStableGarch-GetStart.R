
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
#  FUNCTION:               DESCRIPTION:
#
#  .getStart               Returns initial and boundary values to 
#                          perform optimization 
################################################################################


.getStart <- function(data,m,n,p,q, AR = FALSE, MA = FALSE,
                         cond.dist = c("stable", "gev", "t3", "norm", "std", "sstd", "skstd", "ged"), 
                         TOLG = 1e-7, TOLSTABLE = 2e-2)
{    
  
    # Description:
    #   Get initial values to start the estimation and the bounds for
    #   parameters to be estimated inside GSgarch.Fit function.
    #   Remarks: This function tunes initial parameters to perform optimization
    #   The ARMA coefficients are the ones returned by the "arima" function
    #   adjusted parameters functions from package ("arima" belongs to package "stats" from R)
    #   For GARCH(p,q) Wurtz et al. (2006) suggested
    #   using omega = 0.1, gamma = 0, alpha = 0.1/p and beta = 0.8/q
    #   delta is chosen to be initially equals to 1.5 in almost all the cases 
    #   Keep in mind that delta < alpha for the stable case.
    #   The arma order passed to this function will be 1 even though the 
    #   process does not possess the AR or MA component. Therefore, 
    #   it is not a problem to have AR our MA input parameters equal to FALSE 
    #   even thought the model order says on the contrary.
    
    # Arguments:
    #   data - vector of data
    #   m, n, p, q - model order as in ARMA(m,n)-GARCH/APARCH(p,q)
    #   AR - boolean value that indicates whether we have a model
    #   with th Autoregressive part included
    #   MA - boolean value that indicates whether we have a model
    #   with the Moving Average part included
    #   ARMAonly - Indicates whether we have a pure ARMA model
    #   cond.dist - name of the conditional distribution, one of
    #       gev, stable, norm, std, sstd, ged
    #   TOLSTABLE - boundary tolerance. Should be greater than GSstable.tol
    #   TOLG - pper and lower bounds tolerance. Should be greater than tol
    
    # Return:
    #   Result - A tree columns matrix with columns representing the start, lower and upper
    #     bounds of the parameter set
    
    # FUNCTION:  
    
    # Error control of input parameters
    if (m < 0 || n < 0 || m %% 1 != 0 || n %% 1 != 0)
        stop("'m' and 'n' need to be integers greater than zero")
    
    if (m == 0 || n == 0)
        stop("Expects 'm' and 'n' different from zero")
    
    if ( (m != 1 && AR == TRUE)  || (n != 1 && MA == TRUE) )
        stop("If AR = TRUE, 'm' should be 1 and if MA = TRUE 'n' should be 1")
    
    if (p < 0 || q < 0 || p %% 1 != 0 || q %% 1 != 0)
        stop("'p' and 'q' need to be integers greater than zero")
    
    if(p == 0 && q != 0)
        stop("Invalid Garch(p,q) order")
    
    if( !is.numeric(data) || !is.vector(data))
      stop("data set must be a numerical one dimensional vector")
    
    # Initial variable declaration
    cond.dist = match.arg(cond.dist)
    cond.dist.list = c("stable", "gev", "t3", "norm", "std", "sstd", "skstd", "ged")
    Mean <- mean(data)
    Var <- var(data)
    Dispersion <- mean(abs(x-Mean))
    arima.fit <- c()
    arima.fit.try <- "empty"  
    arima.m <- m
    arima.n <- n
    if(AR == TRUE) # we don't have the AR part
        arima.m <- 0
    if(MA == TRUE) # we don't have the MA part 
        arima.n <- 0
    
    # Try arima fit function to get initial arma parameters
    try(arima.fit.try <- as.vector(arima(data,order = c(arima.m, 0, arima.n))$coef), silent = TRUE)
    if( is.numeric(arima.fit.try) )
    {   
        arima.fit <- arima.fit.try
        if(AR == FALSE && MA == FALSE)
        {
          ar.init <- arima.fit[1:arima.m]
          ma.init <- arima.fit[(arima.m+1):(arima.m+arima.n)]
        }
        if(AR == TRUE && MA == FALSE)
        {
          ar.init <- 0
          ma.init <- arima.fit[(arima.m+1):(arima.m+arima.n)]
        }
        if(AR == FALSE && MA == TRUE)
        {
          ar.init <- arima.fit[1:arima.m]
          ma.init <- 0
        }
        if(AR == TRUE && MA == TRUE)
        {
          ar.init <- 0
          ma.init <- 0
        }
        mean.init <- arima.fit[arima.m+arima.n+1]    
    } else {   
        mean.init <- Mean
        ar.init <- rep(0,m)
        ma.init <- rep(0,n)
    }
 
    # START VALUES
    
    arma.start = c(ar.init, ma.init)
    mu.start = mean.init
    gm.start = rep(0, p)
    
    omega.start = list(
      "stable" = 0.1 * Dispersion,
      "gev" = 0.1 * Var,
      "t3" = 0.1 * Var,
      "norm" = 0.1 * Var,
      "std" = 0.1 * Var,
      "sstd" = 0.1 * Var,
      "skstd" = 0.1 * Var,
      "ged" = 0.1 * Var)
    
    alpha.start = list(
      "stable" = rep(0.1/p, p),
      "gev" = rep(0.1/p, p),
      "t3" = rep(0.1/p, p),
      "norm" = rep(0.1/p, p),
      "std" = rep(0.1/p, p),
      "sstd" = rep(0.1/p, p),
      "skstd" = rep(0.1/p, p),
      "ged" = rep(0.1/p, p))
    
    beta.start = list(
      "stable" = rep(0.8/q, q),
      "gev" = rep(0.4/q, q),
      "t3" = rep(0.8/q, q),
      "norm" = rep(0.8/q, q),
      "std" = rep(0.8/q, q),
      "sstd" = rep(0.8/q, q),
      "skstd" = rep(0.8/q, q),
      "ged" = rep(0.8/q, q))
    
    delta.start = list(
      "stable" = 1.05,
      "gev" = 2,
      "t3" = 2,
      "norm" = 2,
      "std" = 2,
      "sstd" = 2,
      "skstd" = 2,
      "ged" = 2)
    
    skew.start = list(
      "stable" = 0,
      "gev" = 1,
      "t3" = 1,
      "norm" = 1,
      "std" = 1,
      "sstd" = 1,
      "skstd" = 1,
      "ged" = 1)
    
    shape.start = list(
      "stable" = 1.9,
      "gev" = -0.3,
      "t3" = c(2, 4),
      "norm" = 1,
      "std" = 4,
      "sstd" = 4,
      "skstd" = 4,
      "ged" = 4)
        
    # LOWER BOUNDS
    
    mu.lower = - 100 * abs ( mu.start )
    arma.lower = rep ( - 100, m + n )
    omega.lower = TOLG
    alpha.lower = rep ( TOLG, p )
    beta.lower = rep ( TOLG, q )
    gm.lower = rep( - 1 + TOLG, p)
    
    delta.lower = list(
      "stable" = 1,
      "gev" = TOLG,
      "t3" = TOLG,
      "norm" = TOLG,
      "std" = TOLG,
      "sstd" = TOLG,
      "skstd" = TOLG,
      "ged" = TOLG)
    
    skew.lower = list(
      "stable" = - 1 + TOLSTABLE,
      "gev" = 0,
      "t3" = TOLG,
      "norm" = 0,
      "std" = 0,
      "sstd" = TOLG,
      "skstd" = TOLG,
      "ged" = 0)
    
    shape.lower = list(
      "stable" = 1 + TOLSTABLE,
      "gev" = - 100,
      "t3" = c ( TOLG, TOLG),
      "norm" = 0,
      "std" = 2 + TOLG,
      "sstd" = 2 + TOLG,
      "skstd" = 2 + TOLG,
      "ged" = TOLG)
       
    # UPPER BOUNDS
    
    mu.upper = + 100 * abs ( mu.start )
    arma.upper = rep ( + 100, m + n )
    omega.upper = 1000
    alpha.upper = rep ( 1 - TOLG, p )
    beta.upper = rep ( 1 - TOLG, q )
    gm.upper = rep (1 - TOLG, p)

    delta.upper = list(
      "stable" = 2 - TOLG,
      "gev" = 100,
      "t3" = 100,
      "norm" = 100,
      "std" = 100,
      "sstd" = 100,
      "skstd" = 100,
      "ged" = 100)
    
    skew.upper = list(
      "stable" = 1 - TOLSTABLE,
      "gev" = 2,
      "t3" = 100,
      "norm" = 2,
      "std" = 2,
      "sstd" = 100,
      "skstd" = 100,
      "ged" = 2)
    
    shape.upper = list(
      "stable" = 2 - TOLSTABLE,
      "gev" = 100,
      "t3" = c ( 100, 100),
      "norm" = 2,
      "std" = 100,
      "sstd" = 100,
      "skstd" = 100,
      "ged" = 100)    

    # CHECK IF THE STARTING MODELS IS STATIONARY
    start.model.persistency = 1
    start.model.persistency = sum( alpha.start[[cond.dist]] * gsMomentAparch( 
          cond.dist = cond.dist, shape = shape.start[[cond.dist]], 
          skew = skew.start[[cond.dist]], delta = delta.start[[cond.dist]], gm = 0) ) + 
          sum ( beta.start[[cond.dist]] )
    # The real condition is 1 but we let the user to set the limit between 0.95 and 0.999...
    if(start.model.persistency >= 0.95)
    {
      print(start.model.persistency)
      print(paste("The starting model with conditional",cond.dist.list[i], " is not stationary."))
      stop("Change the starting value of the parameters for the reported model.")
    }
    
    # CONSTRUCT THE RESULT
    
    start = c ( mu.start, arma.start, omega.start[[cond.dist]], alpha.start[[cond.dist]],
               gm.start, beta.start[[cond.dist]], delta.start[[cond.dist]], 
               skew.start[[cond.dist]], shape.start[[cond.dist]])
  
    lower = c ( mu.lower, arma.lower, omega.lower, alpha.lower,
              gm.lower, beta.lower, delta.lower[[cond.dist]], 
              skew.lower[[cond.dist]], shape.lower[[cond.dist]])
  
    upper = c ( mu.upper, arma.upper, omega.upper, alpha.upper,
              gm.upper, beta.upper, delta.upper[[cond.dist]], 
              skew.upper[[cond.dist]], shape.upper[[cond.dist]])
    
    namesStart = c("mu", paste("ar", 1:m, sep = ""), 
                   paste("ma", 1:n, sep = ""),
                   "omega", paste("alpha", 1:p, sep = ""), 
                   paste("gm", 1:p, sep = ""), 
                   paste("beta", 1:q, sep = ""), "delta","skew",
                   paste("shape", 1:length(shape.start[[cond.dist]]), sep = ""))        
    
    # Create result
    result = rbind(start,lower,upper)
    colnames(result) = namesStart
    
    # Return
    result
}



# ------------------------------------------------------------------------------



.getStartOld <- function(data,m,n,p,q, AR = FALSE, MA = FALSE, ARMAonly = FALSE,
                         cond.dist = c("stable", "gev", "t3", "norm", "std", "sstd", "skstd", "ged"), 
                         TOLG = 1e-7, TOLSTABLE = 2e-2)
{    
  # Description:
  #   Get initial values to start the estimation and the bounds for
  #   parameters to be estimated inside GSgarch.Fit function.
  #   Remarks: This function tunes initial parameters to perform optimization
  #   The ARMA coefficients are the ones returned by the "arima" function
  #   adjusted parameters functions from package ("arima" belongs to package "stats" from R)
  #   For GARCH(p,q) Wurtz et al. (2006) suggested
  #   using omega = 0.1, gamma = 0, alpha = 0.1/p and beta = 0.8/q
  #   delta is chosen to be initially equals to 1.5 in almost all the cases 
  #   Keep in mind that delta < alpha for the stable case.
  #   The arma order passed to this function will be 1 even though the 
  #   process does not possess the AR or MA component. Therefore, 
  #   it is not a problem to have AR our MA input parameters equal to FALSE 
  #   even thought the model order says on the contrary.
  
  # Arguments:
  #   data - vector of data
  #   m, n, p, q - model order as in ARMA(m,n)-GARCH/APARCH(p,q)
  #   AR - boolean value that indicates whether we have a model
  #   with th Autoregressive part included
  #   MA - boolean value that indicates whether we have a model
  #   with the Moving Average part included
  #   ARMAonly - Indicates whether we have a pure ARMA model
  #   cond.dist - name of the conditional distribution, one of
  #       gev, stable, norm, std, sstd, ged
  #   TOLSTABLE - boundary tolerance. Should be greater than GSstable.tol
  #   TOLG - pper and lower bounds tolerance. Should be greater than tol
  
  # Return:
  #   Asdf - The asdf 
  
  # FUNCTION:  
  
  # Error control
  if (m < 0 || n < 0 || m %% 1 != 0 || n %% 1 != 0)
    stop("'m' and 'n' need to be integers greater than zero")
  if (m == 0 || n == 0)
    stop("GSgarch.GetStart expects 'm' and 'n' different from zero")
  if ( (m != 1 && AR == TRUE)  || (n != 1 && MA == TRUE) )
    stop("If AR = TRUE, 'm' should be 1 and if MA = TRUE 'n' should be 1")
  if (p < 0 || q < 0 || p %% 1 != 0 || q %% 1 != 0)
    stop("'p' and 'q' need to be integers greater than zero")
  if(p == 0 && q != 0)
    stop("Invalid Garch(p,q) order")
  # Conditional distribution
  cond.dist = match.arg(cond.dist)
  if( !is.numeric(data) || !is.vector(data))
    stop("data set must be a numerical one dimensional vector")
  
  # Initial ARMA parameters
  Mean <- mean(data)
  Var <- var(data)
  Dispersion <- mean(abs(x-Mean))
  arima.fit <- c()
  arima.failed <- FALSE
  arima.fit.try <- "empty"  
  arima.m <- m
  arima.n <- n
  if(AR == TRUE) # we don't have the AR part
    arima.m <- 0
  if(MA == TRUE) # we don't have the MA part 
    arima.n <- 0
  
  if(ARMAonly == FALSE) # try 'arima' function only if we are not in a pure 'arma' model
  {
    # try arima fit with correct order depending on AR and MA
    try(arima.fit.try <- as.vector(arima(data,order = c(arima.m, 0, arima.n))$coef), silent = TRUE)
    if( is.numeric(arima.fit.try) )
      arima.fit <- arima.fit.try
    else
      arima.failed <- TRUE
    
    # Configuring initial ar,ma and mean parameters returned by 'arima' function.
    if (arima.failed == FALSE) 
    {
      if(AR == FALSE && MA == FALSE)
      {
        ar.init <- arima.fit[1:arima.m]
        ma.init <- arima.fit[(arima.m+1):(arima.m+arima.n)]
      }
      if(AR == TRUE && MA == FALSE)
      {
        ar.init <- 0
        ma.init <- arima.fit[(arima.m+1):(arima.m+arima.n)]
      }
      if(AR == FALSE && MA == TRUE)
      {
        ar.init <- arima.fit[1:arima.m]
        ma.init <- 0
      }
      if(AR == TRUE && MA == TRUE)
      {
        ar.init <- 0
        ma.init <- 0
      }
      mean.init <- arima.fit[arima.m+arima.n+1]
    }
    else # arima function failed
    {
      mean.init <- Mean
      ar.init <- rep(0,m)
      ma.init <- rep(0,n)
    }
  } else {
    mean.init <- Mean
    ar.init <- rep(0,m)
    ma.init <- rep(0,n)        
  }
  
  # Initial APARCH and Density parameters
  omega.init <- min(0.1,0.1*Var)
  alpha.init <- rep(0.1/p,p)
  beta.init <- rep(0.8/q,q)
  shape.init <- 1
  skew.init <- 1
  gm.init <- rep(0,p)
  delta.init <- 1.5
  sigma.init <- Var
  
  # Lower Bound
  #mean.lower <- mean.init - 3*abs(mean.init)
  mean.lower <- -10*abs(mean.init)
  arma.lower <- rep(-10+TOLG,m+n)
  omega.lower <- rep(TOLG,1)
  alpha.lower <- rep(TOLG,p)
  beta.lower <- rep(TOLG,q)
  shape.lower <- 0
  skew.lower <- 0
  gm.lower <- rep(-1 + TOLG,p)
  delta.lower <- TOLG
  sigma.lower <- 0
  
  # Upper Bound  
  #mean.upper <- mean.init + 3*abs(mean.init)
  mean.upper <- 10*abs(mean.init)
  arma.upper <- rep(10-TOLG,m+n)
  omega.upper <- rep(1-TOLG,1)
  #omega.upper <- 100*var(data)
  alpha.upper <- rep(1-TOLG,p)
  beta.upper <- rep(1-TOLG,q)
  shape.upper <- 4
  skew.upper <- 4
  gm.upper <- rep(1 - TOLG,p)
  # delta.upper <- 3 + TOLG
  delta.upper <- 30 + TOLG
  sigma.upper <- 10*Var
  
  # Setting skew and shape appropriatelly for other conditional Distributions
  if (cond.dist == "std")
  {   
    shape.init <- 4
    shape.lower <- 2 + TOLG; shape.upper <- 20   
  }
  if (cond.dist == "sstd")
  {   
    shape.init <- 4; skew.init = 1
    shape.lower <- 2 + TOLG; shape.upper <- 20
    skew.lower <- TOLG; skew.upper <- 20 
  }
  if (cond.dist == "skstd")
  {   
    shape.init <- 4; skew.init = 1
    shape.lower <- TOLG; shape.upper <- 100
    skew.lower <- TOLG; skew.upper <- 100   
  }
  if (cond.dist == "ged")
  {   
    shape.init <- 4
    shape.lower <- 0 + TOLG; shape.upper <- 20   
  }
  if (cond.dist == "t3")
  {   
    shape.init <- c(4,1)
    shape.lower <- rep(0 + TOLG,2)
    shape.upper <- rep(100,2)   
  }
  if(cond.dist == "gev")
  {
    mean.init <- 0
    shape.init <- 0.1; 
    alpha.init <- rep(0.05/p,p)
    beta.init <- rep(0.8/q,q)
    # shape.lower <- -0.5 + TOLG; shape.upper <- 0.5 - TOLG
    shape.lower <- -10 + TOLG; shape.upper <- 10 - TOLG
  }
  if(cond.dist == "stable")
  {
    shape.init <- 1.8; delta.init <- 1
    shape.lower <- 1 + TOLSTABLE; shape.upper <- 2 - TOLSTABLE
    skew.init <- 0
    skew.lower <- -1 + TOLSTABLE; skew.upper <- 1 - TOLSTABLE
    delta.upper <- 1.9
  }
  
  # Constructing the initial, lower and upper bound vectors
  arma.init <- c(ar.init,ma.init)
  if(!ARMAonly)
  {
    init <- c(mean.init,arma.init,omega.init,alpha.init,gm.init,
              beta.init,delta.init,skew.init,shape.init)
    lower <- c(mean.lower,arma.lower,omega.lower,alpha.lower,gm.lower,
               beta.lower,delta.lower,skew.lower,shape.lower)
    upper <- c(mean.upper,arma.upper,omega.upper,alpha.upper,gm.upper,
               beta.upper,delta.upper,skew.upper,shape.upper)
    
    namesStart = c("mu", paste("ar", 1:length(ar.init), sep = ""), paste("ma", 1:length(ma.init), sep = ""),
                   "omega", paste("alpha", 1:length(alpha.init), sep = ""),paste("gm", 1:length(gm.init), sep = ""), 
                   paste("beta", 1:length(beta.init), sep = ""), "delta","skew",paste("shape", 1:length(shape.init), sep = ""))        
    
  } else {
    init <- c(mean.init,arma.init,skew.init,shape.init,sigma.init)
    lower <- c(mean.lower,arma.lower,skew.lower,shape.lower,sigma.lower)
    upper <- c(mean.upper,arma.upper,skew.upper,shape.upper,sigma.upper)   
    namesStart = c("mu", paste("ar", 1:length(ar.init), sep = ""), paste("ma", 1:length(ma.init), sep = ""),
                   "skew",paste("shape", 1:length(shape.init), sep = ""),"sigma")
  }
  
  if(arima.failed == TRUE)
    warning("arima function from package stats failed to get initial AR and MA coefficients")
  
  # Create result
  result = rbind(init,lower,upper)
  colnames(result) = namesStart
  
  # Return
  result
}







################################################################################
