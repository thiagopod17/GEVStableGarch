
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
#  GSgarch.GetStart        Returns initial and boundary values to 
#                          perform optimization 
################################################################################


GSgarch.GetStart <- function(data,m,n,p,q, AR = FALSE, MA = FALSE,ARMAonly = FALSE,
                             cond.dist = "norm", GSstable.tol.b = 2e-2, GStol.b = 1e-7)
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
    #   GSstable.tol.b - boundary tolerance. Should be greater than GSstable.tol
    #   GStol.b - pper and lower bounds tolerance. Should be greater than tol

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
    cond.dist.list <- c("norm", "std", "sstd", "gev", "stable","ged")
    if( !any(cond.dist.list == cond.dist) )   
        stop ("Invalid Conditional Distribution. Choose: norm,std,sstd,gev,ged or stable")
    if( !is.numeric(data) || !is.vector(data))
        stop("data set must be a numerical one dimensional vector")
    
    # Initial ARMA parameters
    Mean <- mean(data)
    Var <- var(data)
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
    arma.lower <- rep(-10+GStol.b,m+n)
    omega.lower <- rep(GStol.b,1)
    alpha.lower <- rep(GStol.b,p)
    beta.lower <- rep(GStol.b,q)
    shape.lower <- 0
    skew.lower <- 0
    gm.lower <- rep(-1 + GStol.b,p)
    delta.lower <- GStol.b
    sigma.lower <- 0
    
    # Upper Bound	
    #mean.upper <- mean.init + 3*abs(mean.init)
    mean.upper <- 10*abs(mean.init)
    arma.upper <- rep(10-GStol.b,m+n)
    omega.upper <- rep(1-GStol.b,1)
    #omega.upper <- 100*var(data)
    alpha.upper <- rep(1-GStol.b,p)
    beta.upper <- rep(1-GStol.b,q)
    shape.upper <- 4
    skew.upper <- 4
    gm.upper <- rep(1 - GStol.b,p)
    delta.upper <- 3 + GStol.b
    sigma.upper <- 10*Var
    
    # Setting skew and shape for other conditional Distributions
    if (cond.dist == "std")
    {   
        shape.init <- 4
        shape.lower <- 2 + GStol.b; shape.upper <- 20   
    }
    if (cond.dist == "sstd")
    {   
        shape.init <- 4
        shape.lower <- 2 + GStol.b; shape.upper <- 20
        skew.lower <- GStol.b; skew.upper <- 20   
    }
    if (cond.dist == "ged")
    {   
      shape.init <- 4
      shape.lower <- 0 + GStol.b; shape.upper <- 20   
    }
    if(cond.dist == "gev")
    {
        mean.init <- 0
        shape.init <- 0.01; beta.init <- rep(0.4/q,q)
        shape.lower <- -0.5 + GStol.b; shape.upper <- 0.5 - GStol.b
    }
    if(cond.dist == "stable")
    {
        shape.init <- 1.8; delta.init <- 1
        shape.lower <- 1 + GSstable.tol.b; shape.upper <- 2 - GSstable.tol.b
        skew.init <- 0
        skew.lower <- -1 + GSstable.tol.b; skew.upper <- 1 - GSstable.tol.b
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
    } else {
      init <- c(mean.init,arma.init,skew.init,shape.init,sigma.init)
      lower <- c(mean.lower,arma.lower,skew.lower,shape.lower,sigma.lower)
      upper <- c(mean.upper,arma.upper,skew.upper,shape.upper,sigma.upper)      
    }
      
    if(arima.failed == TRUE)
        warning("arima function from package stats failed to get initial AR and MA coefficients")
    
    # return argument
    return(rbind(init,lower,upper))
}



################################################################################
