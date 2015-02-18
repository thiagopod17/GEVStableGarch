
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
	  if (spec@distribution == "sstd") 
	  z = rsstd(n, nu = model$shape, xi = model$skew)

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
    if(!is.null(spec@model$alpha)){ # Arch, Garch or Aparch model
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
	} else { # pure ARMA, AR or MA model.
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