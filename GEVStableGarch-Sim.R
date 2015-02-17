
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
# FUNCTION:               SIMULATION:
#  GSgarch.Sim            Simulates a GARCH/APARCH process with GEV or stable
#						  conditional distribution
################################################################################
library(stabledist)


GSgarch.Sim <-
    function(spec = garchSpec(), n = 100, n.start = 100)
{
    # Description:
    #   Simulates a time series process from the GARCH family

    # Arguments:
    #   model - a specification object of class 'fGARCHSPEC' as
    #     returned by the function \code{garchSpec}:
    #     ar - a vector of autoregressive coefficients of
    #       length m for the ARMA specification,
    #     ma - a vector of moving average coefficients of
    #       length n for the ARMA specification,
    #     omega - the variance value for GARCH/APARCH
    #       specification,
    #     alpha - a vector of autoregressive coefficients
    #       of length p for the GARCH/APARCH specification,
    #     gamma - a vector of leverage coefficients of
    #       length p for the APARCH specification,
    #     beta - a vector of moving average coefficients of
    #       length q for the GARCH/APARCH specification,
    #     mu - the intercept for ARMA specification (mean=mu/(1-sum(ar))),
    #     delta - the exponent value used in the variance
    #       equation.
    #     skew - a numeric value for the skew parameter.
    #     shape - a numeric value for the shape parameter.
    #   n - an integer, the length of the series
    #   n.start - the length of the warm-up sequence to reduce the
    #       effect of initial conditions.

    # FUNCTION:

	# Error treatment of input parameters
	if(n < 2)
	   stop("The parameter 'n' must be > 2")

    # Specification:
    stopifnot(class(spec) == "fGEVSTABLEGARCHSPEC")
    model = spec@model

    # Random Seed:
    if (spec@rseed != 0) set.seed(spec@rseed)

    # Enlarge Series:
    n = n + n.start

    # Create Innovations:
    if (spec@distribution == "gev")
        z <- rgev(n, xi = model$shape)
    if (spec@distribution == "stable")
        z <- stabledist::rstable(n = n, alpha = model$shape, beta = model$skew, pm = 2)
    if (spec@distribution == "norm")
        z = rnorm(n)
    if (spec@distribution == "std")
        z = rstd(n, nu = model$shape)    

    # Expand to whole Sample:    NAO ENTENDI PORQUE USAR A FUNCAO rev()???
    delta = model$delta
    if(!is.null(spec@model$alpha)){
    	z = c(rev(spec@presample[, 1]), z)
    	h = c(rev(spec@presample[, 2]), rep(NA, times = n))
    	y = c(rev(spec@presample[, 3]), rep(NA, times = n))
    	m = length(spec@presample[, 1])
    	names(z) = names(h) = names(y) = NULL
    } else {
    	z = c(rev(spec@presample[, 1]), z)
    	y = c(rev(spec@presample[, 2]), rep(NA, times = n))
    	m = length(spec@presample[, 1])
    	names(z) = names(y) = NULL 	
    }
    
    # Determine Coefficients:
    mu = model$mu
    ar = model$ar
    ma = model$ma
    omega = model$omega
    alpha = model$alpha
    gamma = model$gamma
    beta = model$beta
    deltainv = 1/delta

    # Determine Orders:
    order.ar = length(ar)
    order.ma = length(ma)
    order.alpha = length(alpha)
    order.beta = length(beta)

    # Iterate GARCH / APARCH Model and create Sample:
    if(!is.null(spec@model$alpha)){
    	eps = h^deltainv*z   # here the variable 'h' represents the process '(sigma_t)^delta'
    	for (i in (m+1):(n+m)) {
       	 	h[i] =  omega +
            	sum(alpha*(abs(eps[i-(1:order.alpha)]) -
                gamma*(eps[i-(1:order.alpha)]))^delta) +
            	sum(beta*h[i-(1:order.beta)])
        	eps[i] = h[i]^deltainv * z[i]
        	y[i] = mu  +
            	sum(ar*y[i-(1:order.ar)]) +
            	sum(ma*eps[i-(1:order.ma)]) + eps[i]
    	}
    	# Sample:
    	data = cbind(
       	 	z = z[(m+1):(n+m)],
        	sigma = h[(m+1):(n+m)]^deltainv,
        	y = y[(m+1):(n+m)])    	
	} else {
		eps = z
		for (i in (m+1):(n+m)) {
        	y[i] = mu  +
            	sum(ar*y[i-(1:order.ar)]) +
            	sum(ma*eps[i-(1:order.ma)]) + eps[i]
    	}
    	data = cbind(
       	 	z = z[(m+1):(n+m)],
        	y = y[(m+1):(n+m)])  
	}
    
    rownames(data) = as.character(1:n)
    if(n.start > 0)
    	data = data[-(1:n.start),]


    # Return Values:
    from <-
        timeDate(format(Sys.time(), format = "%Y-%m-%d")) - NROW(data)*24*3600
    charvec  <- timeSequence(from = from, length.out = NROW(data))
    if(!is.null(spec@model$alpha)){
    	ans <- timeSeries(data = data[, c(3,2,1)], charvec = charvec)
    	colnames(ans) <- c("garch", "sigma", "eps")
	} else {
		ans <- timeSeries(data = data[, c(2,1)], charvec = charvec)		
	    colnames(ans) <- c("series", "eps")
	}
    
    attr(ans, "control") <- list(garchSpec = spec)

    # Return Value:
    ans
}


################################################################################



################################################################################
# TEST CASES FOF FUNCTION:               SPECIFICATION:
#  GSgarch.Sim               Useful for finding erros. Please add the expected 
#							 output of the test case whenever possible. 
#							 This will help to trace erros during debugging.
################################################################################
library("fGarch")
library("fExtremes")
library("stabledist")
library("skewt")
library("Rsolnp")

# Simulate AR(1)-GARCH(1,1) with conditional GEV distr.
spec <- GSgarchSpec(model = list(ar = 0.1, alpha = 0.3, beta = 0.2, delta = 2), 
presample = NULL,cond.dist = "gev",rseed = NULL)  
sim <- GSgarch.Sim(spec, n = 1000, n.start = 0,
             extended = TRUE)
plot(sim)
# Simulate GARCH(1,1) with conditional stable distr.
spec <- GSgarchSpec(model = list(alpha = 0.05, beta = 0.01, omega = 0.01, delta = 2,shape = 1.2), 
presample = NULL,cond.dist = "stable",rseed = NULL)  
sim <- GSgarch.Sim(spec, n = 100, n.start = 0)
plot(sim)
# Simulate with small sample size: Expect error about size of n.
spec <- GSgarchSpec(model = list(alpha = 0.05, beta = 0.01, omega = 0.01, delta = 2,shape = 1.2), 
presample = NULL,cond.dist = "stable",rseed = NULL)  
sim <- GSgarch.Sim(spec, n = 1, n.start = 0)
# simulate pure ARMA(1,1) with conditional stable distribution
spec <- GSgarchSpec(model = list(ar = 0.05, ma = 0.01,shape = 1.5,skew = -0.5), 
presample = NULL,cond.dist = "stable",rseed = NULL)  
sim <- GSgarch.Sim(spec, n = 100, n.start = 0)
plot(sim)
# simulate pure ARMA(2,2) with conditional stable distribution
spec <- GSgarchSpec(model = list(ar = c(0.05,0.3), ma = c(0.01,0.02),shape = 0.5,skew = -0.5), 
presample = NULL,cond.dist = "stable",rseed = NULL)  
sim <- GSgarch.Sim(spec, n = 100, n.start = 0)
plot(sim)























############
# Simulation 
############
GSgarch.Sim <- function(N = 1000,mu = 0.1,a = c(0.5,0.3),
b = c(-0.4,0.3,-0.1), omega = 0.05, alpha = c(0.1), gm = c(0),
beta = c(0.1,0.05,0.03), delta = 2, skew = 0, shape = 3, cond.dist = "norm")
{
	# only alpha different from zero.
	
	
	
	# error treatment on input parameters
    cond.dist.list <- c("norm", "std", "sstd", "gev", "stable")
    if( !any(cond.dist.list == cond.dist) )   
        stop ("Invalid Distribution. 
        Choose: norm,std,sstd,gev or stable")
    if(length(alpha) == 0)
        stop("'alpha' must be a non zero vector")
    if(length(delta) != 0 ||  length(gm) != 0) # means aparch model
    {
      if(length(alpha) != length(gm) || length(delta) != 1)
        stop("'alpha' and 'gm' must have the same size for APARCH models")
    }
    if(sum(!(abs(gm)<1)) > 0)
        stop("'gm[i]' should satisfy -1 < gm < 1")
        
        
	# variable declaration
    x <- rep(mu,N); e <- rep(0,N); h <- rep(0.1,N); z <- rep(0,N)
    m <- length(a); n <- length(b); mn <- max(m,n)
    p <- length(alpha); q <- length(beta); pq <- max(p,q); 
    mnpq <- max(mn,pq); lambda <- 1;
	  if(!(N > 0) || (N < 1+mnpq+1))
	    stop("invalid size simulate. 'N' should satisfy N > 2+max(m,n,p,q)")

    # Conditional distribution
    if (cond.dist == "norm")
        z <- rnorm(N)
    if (cond.dist == "std")
        z <- rstd(N, nu = shape)
    if (cond.dist == "gev")
        z <- rgev(N, xi = shape)
    if (cond.dist == "stable")
    {
        z <- stabledist::rstable(N, alpha = shape, beta = skew, pm = 2)
        lambda <- mean((abs(z) + gm*z)^delta)
    }
    # we will not verify stationarity conditions to allow any time series to be simulated
    x[1:mn] <- 0; e[1:pq] <- 0
    for (t in (1+mnpq):(N)){
        h[t] <- omega + 
                sum(alpha*(abs(e[(t-1):(t-p)])-gm*e[(t-1):(t-p)])^delta) + 
                sum(beta*h[(t-1):(t-q)])
        e[t] <- z[t]*(h[t])^(1/delta)
        x[t] <- mu + sum(a*x[(t-1):(t-m)]) + 
                sum(b*e[(t-1):(t-n)]) + e[t]
    }
    garch <- matrix(NA, nrow = N, ncol = 2)
    garch[,1] <- x
    garch[,2] <- h^(1/delta)
	 dimnames(garch) <- list(rownames(garch, do.NULL = FALSE, prefix = ""),
	                     c("x[t]","sigma[t]"))
	 out.sim <- NULL # output of the GSgarch.Sim function
	 TMPvector <- c(if(m != 0 || n != 0) "arma(",m,",",n,")-",
	                 if(length(gm)==0)c("garch(",p,",",q,")"),
	                 if(length(gm)>0)c("aparch(",p,",",q,")"))
	 TMPorder <- paste(TMPvector, collapse="")
	 TMPvectorintercept <- c("Intercept:",if(length(mu)==0)"TRUE", 
	                    if(length(mu)==1)"FALSE")
	 TMPintercept <- paste(TMPvectorintercept, collapse="")
	 out.sim$model <- paste(TMPorder,"##",TMPintercept, collapse="")
	 out.sim$cond.dist <- cond.dist
   out.sim$series <- garch
return(out.sim)
}
