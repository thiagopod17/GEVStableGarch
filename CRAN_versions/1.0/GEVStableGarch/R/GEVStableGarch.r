###################################################################
# Package Name: GEVStableGarch 
# Authors: Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano
# and Silvia Regina Costa Lopes
# Date: 2014/03/08
# License: The functions may be distributed free of charge and 
# used by anyone if credit is given. It has been tested, 
# but it comes with no guarantees and the authors assume no 
# liability for its use or misuse.
###################################################################


###########################
###  Time series model  ###
# According to Wurtz et al. (2006) 
# xt = µ + a(B)xt + b(B)et,
# et = zt * sigmat
# zt ~ Dv(0, 1)
# sigmat^2 = omega + sum(alphai(e(t-i))^2) + sum(betaj*sigmat-j^2)

#################################################
###  Variables notation used inside functions ###
# N: sample size
# m,n: ARMA order
# p,q, pq = max(p,q): GARCH order. They were called u,v, uv in Wurtz et al.(2006)
# mu, a, b: ARMA parameters
# omega, alpha, gamma, beta: APARCH vector parameters
# garchLLH: Log Likelihood for the ARMA(m,n)-GARCH(p,q) model
# parm: ARMA-GARCH parameters concatenated as c(mu,a,b,omega,alpha,gamma,beta)
# h: conditional variance. It was referenced as sigmat in model equations
# x: data set. It is the time series that will be modelled as ARMA-GARCH
# xi: GEV shape parameter. xi > -0.5 (llh) and < 0.5 (finiteness of variance, see Zhao et al. (2011))
# AR, MA: if AR (or MA) equals TRUE, then the model has AR (or MA) order equals to zero.
# param or pm (stabledist): Parametrization of stable distributions as chosen to be 2.  
# Even for stable distribution llh, there's a problem for finding the estimators for
# alpha near to 2 and beta > 0. The bad performance on the ARMA-GARCH model was
# similar to the single llh estimation of i.i.d samples of stable distribution. 
# In our simulations, the param = 2 was chosen to be the ones that performs
# better in these situations. 

######################################
### Constants to adjust tolerances ###
# GStol <- 1e-8 # general tolerance for arma-garch parameters. In the beggining it was set to 1e-5
# GStol.b <- 1e-7 # upper and lower bounds tolerance. Should be greater than tol
# GSstable.tol <- 1e-2 # tolerance for stable distribution parameter set
# GSstable.tol.b <- 2e-2 # boundary tolerance. Should be greater than GSstable.tol

#################
#### Remarks ####
# et for ARMA(m,n), n > 1 will be initiated as 0, i.e, e[-n+1:0] = 0.
# ht will be initiated as 0.1 (see eq. (22) of Wurtz, 2006).

######################################################################################
# FUNCTIONS IN THIS FILE
######################################################################################
#  NAME   	 		  DESCRIPTION
#  GSgarch.GetStart	  Returns initial and boundary values to perform optimization 
#  GSgarch.Dist		  Computes density values for several conditional distributions
#  GSgarch.Fit		  Fits ARMA-GARCH or ARMA-APARCH models
#  GSgarch.Sim		  Simulate ARMA-GARCH or ARMA-APARCH models
#  GSgarch.GetOrder   Return a matrix with parameter order for use in GSgarch.FitAIC
#  GSgarch.FitAIC	  Find best fitted model according to AIC criterion
######################################################################################


########################################################
# Initial configurations for stable distribution
# computation
########################################################


.onAttach <- function(libname, pkgname)
{
    options(.stableIsLoaded=FALSE)
	if(length(find.package("stable",quiet = TRUE)))
		options(.stableIsLoaded=TRUE)	
    if(getOption('.stableIsLoaded', default = FALSE) == TRUE)
    {
        GSgarch.dstable <<- function(x,alpha,beta = 0, gamma = 1, 
        delta = 0, param = 0)
	  {
	       return(stable::dstable.quick(x, alpha, beta, gamma, 
             delta, param))
	  }
    }
    if(getOption('.stableIsLoaded', default = FALSE) == FALSE)
    {
        GSgarch.dstable <<- function(x,alpha,beta = 0, gamma = 1, 
        delta = 0, param = 0)
	  {
	       return(stabledist::dstable(x, alpha, beta, gamma, 
		 delta, pm = param))
	  }
    }
}
########################################################
# get initial, lower and upper bound for the parameters 
# to perform optimization
########################################################
GSgarch.GetStart <- function(data,m,n,p,q, AR = FALSE, MA = FALSE, cond.dist = "norm", 
					GSstable.tol.b = 2e-2, GStol.b = 1e-7)
{
# Remarks: This function tunes initial parameters to perform optimization
# The ARMA coefficients are the ones returned by the "arima" function
# adjusted parameters functions from package ("arima" belongs to package "stats" from R)
# For GARCH(p,q) Wurtz et al. (2006) suggested
# using omega = 0.1, gamma = 0, alpha = 0.1/p and beta = 0.8/q
# delta is chosen to be initially equals to 1.5 in almost all the cases 
# Keep in mind that delta < alpha for the stable case.
# The arma order passed to this function will be 1 even though the 
# process does not possess the AR or MA component. Therefore, 
# it is not a problem to have AR our MA input parameters equal to FALSE 
# even thought the model order says on the contrary.
###
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
	cond.dist.list <- c("norm", "std", "sstd", "gev", "stable")
  if( !any(cond.dist.list == cond.dist) )   
    stop ("Invalid Conditional Distribution. Choose: norm,std,sstd,gev or stable")
	if( !is.numeric(data) || !is.vector(data))
		stop("data set must be a numerical one dimensional vector")
    # Initial ARMA parameters
    arima.fit <- c()
    arima.failed <- FALSE
    arima.fit.try <- "empty"  
    arima.m <- m
    arima.n <- n
    if(AR == TRUE) # we don't have the AR part
	  arima.m <- 0
    if(MA == TRUE) # we don't have the MA part 
	  arima.n <- 0
    # try arima fit with correct order depending on AR and MA
	try(arima.fit.try <- as.vector(arima(data,order = c(arima.m, 0, arima.n))$coef), silent = TRUE)
	if( is.numeric(arima.fit.try) )
	  arima.fit <- arima.fit.try
    else
	  arima.failed <- TRUE
    if (arima.failed == FALSE) # we have initial parameters returned by 'arima' function.
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
    mean.init <- mean(data)
    ar.init <- rep(0,m)
    ma.init <- rep(0,n)
  }
  
	# Initial APARCH and Density parameters
	omega.init <- 0.1
	alpha.init <- rep(0.1/p,p)
	beta.init <- rep(0.8/q,q)
	shape.init <- 1
	skew.init <- 1
	gm.init <- rep(0,p)
	delta.init <- 1.5

	# Lower Bound
	mean.lower <- mean.init - 3*abs(mean.init)
	arma.lower <- rep(-10+GStol.b,m+n)
	omega.lower <- rep(GStol.b,1)
	alpha.lower <- rep(GStol.b,p)
	beta.lower <- rep(GStol.b,q)
	shape.lower <- 0
	skew.lower <- 0
	gm.lower <- rep(-1 + GStol.b,p)
	delta.lower <- GStol.b

	# Upper Bound	
	mean.upper <- mean.init + 3*abs(mean.init)
	arma.upper <- rep(10-GStol.b,m+n)
	omega.upper <- rep(1-GStol.b,1)
	alpha.upper <- rep(1-GStol.b,p)
	beta.upper <- rep(1-GStol.b,q)
	shape.upper <- 4
	skew.upper <- 4
	gm.upper <- rep(1 - GStol.b,p)
	delta.upper <- 3 + GStol.b

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
	arma.init <- c(ar.init,ma.init)
	init <- c(mean.init,arma.init,omega.init,alpha.init,gm.init,
			  beta.init,delta.init,skew.init,shape.init)
	lower <- c(mean.lower,arma.lower,omega.lower,alpha.lower,gm.lower,
			   beta.lower,delta.lower,skew.lower,shape.lower)
	upper <- c(mean.upper,arma.upper,omega.upper,alpha.upper,gm.upper,
			   beta.upper,delta.upper,skew.upper,shape.upper)
	if(arima.failed == TRUE)
		warning("arima function from package stats failed to get initial AR and MA coefficients")
    return(rbind(init,lower,upper))
}
##########################
# conditional distribution
##########################
GSgarch.Dist = function(z, hh, shape = 4, skew = 0.1, cond.dist = "sstd", GStol = 1e-8) 
{
    # We choose not to avoid calling this function with out of bound parameters
	  # because in these cases we simply return a huge likelihood to guide the optimization algorithm.
	  cond.dist.list <- c("norm", "std", "sstd", "gev", "stable")
    if( !any(cond.dist.list == cond.dist) )   
        stop ("Invalid Conditional Distribution. Choose: norm,std,sstd,gev or stable")
    if(sum(is.na(hh)) > 0 || min(hh) == 0)
    {
	  #stop ("NA or zero element found in vector hh")
        warning("NA or zero element found in vector hh")
	  return(1e99)
    }

	  # Processing
    if(cond.dist == "norm")
        return(-sum(log(dnorm(x = z/hh)/hh)))
    if(cond.dist == "std")
    {
		  if(!(shape > 2))
			stop("Invalid shape in std distribution. shape > 2")
      nu = shape
      return(-sum(log(dstd(x = z/hh, nu = nu)/hh)))
    }
    if(cond.dist == "sstd")
    {
		  if(!(shape > 2) || !(skew > 0))
			  stop("Invalid shape or skew in skewt parameters. shape > 2 and skew > 0")
      gm = skew
      xi = gm
      nu = shape
	    return(-sum(log(dsstd(x = z/hh, nu = nu, xi = xi)/hh)))
    }
    if(cond.dist == "gev")
    {
		  if(abs(shape) < 1e-6)
			  stop("shape parameter from GEV is to small. abs(shape) < 1e-6")
      sig <- hh
      xi <- shape
      arg <- z/sig
      y <- 1 + xi * arg
      if(sum(is.na(y)) || sum(is.nan(y)) ||
      sum(is.infinite(y))) 
          return(1e99)
      gev.cond.gD <- any(y < GStol ) || abs(xi) < 1e-6
      if(gev.cond.gD)
      {
          return(1e99)
      }
      llh <- sum(log(sig)) + sum(y^(-1/xi)) + sum(log(y))*(1/xi + 1)
      return(llh)
    }
    if(cond.dist == "stable")
    {
		    if( !(shape > 1) || !(shape < 2) || !(abs(skew) < 1))
			      stop("Invalid shape or skew in stable parameters. 1 < shape < 2 and |skew| < 1")
		    sig <- hh
        alpha.stable <- shape
        beta.stable <- skew
        arg <- z/sig
        y <- arg
        if(sum(is.na(y)) || sum(is.nan(y)) ||
        sum(is.infinite(y)))
            return(1e99)
        #dens.stable <- stable::dstable.quick(y,alpha = alpha.stable, 
        #beta = beta.stable, gamma = 1, delta = 0, param = 2)
		dens.stable <- GSgarch.dstable(y,alpha = alpha.stable,
		beta = beta.stable, gamma = 1,delta = 0, param = 2)
        llh <- sum(log(sig[sig>0])) - 
               sum(log(dens.stable[dens.stable>0]))
        return(llh)
    }
}
##################
# ARMA-GARCH fit
##################
GSgarch.Fit <- function(data, m,n,p,q, intercept = TRUE, printRes = FALSE, 
                cond.dist = "norm", APARCH = FALSE, algorithm = "sqp",
                get.res = FALSE, control = NULL, GSstable.tol = 1e-2, GStol = 1e-8)
{  
    # Error Control: Stop if some conditions are not met
    #if (!is.numeric(data) || is.NA(data) || is.NULL(data) || is.Inf(data))
    #    stop("Invalid 'data' input. It may be contain NA, NULL or Inf.")
    cond.dist.list <- c("norm", "std", "sstd", "gev", "stable")
    algoritm.list <- c("sqp","sqp.rest","nlminb")
    if( !any(cond.dist.list == cond.dist) )   
        stop ("Invalid Distribution. Choose: norm,std,sstd,gev or stable")
    if( !any(algoritm.list == algorithm) )   
        stop ("Invalid Algorithm. Choose: sqp, sqp.rest or nlminb")
    if(m%%1 != 0 || n%%1 != 0 || p%%1 != 0 || q%%1 != 0 || 
       any (c(m,n,p,q) < 0) || (p == 0 && q != 0) || (p == 0 && APARCH) ||
       any (c(m,n,p,q) > 10) ) 
       stop ("Invalid ARMA-GARCH order. We allow pure GARCH or APARCH. AR/MA/ARMA-GARCH/APARCH models.
	          The order of the parameters could be set up to 10.")

    data <- data; N <- length(data)
    out <- NULL # output of the GSgarch.Fit function
    out$order <- c(m,n,p,q,intercept,APARCH)
    TMPvector <- c(if(m != 0 || n != 0) c("arma(",m,",",n,")-"),
                   if(APARCH==FALSE)c("garch(",p,",",q,")"),
                   if(APARCH==TRUE)c("aparch(",p,",",q,")"))
    TMPorder <- paste(TMPvector, collapse="")
    TMPvectorintercept <- c("Intercept:",if(intercept==TRUE)"TRUE", 
                            if(intercept==FALSE)"FALSE")
    TMPintercept <- paste(TMPvectorintercept, collapse="")
    out$model <- paste(TMPorder,"##",TMPintercept, collapse="")
    out$cond.dist <- cond.dist
    out$data <- data
    ARMAonly <- FALSE
    optim.finished <- FALSE
    AR = FALSE; MA <- FALSE; GARCH <- FALSE
    if( m == 0) AR <- TRUE
    if( n == 0) MA <- TRUE
    if( q == 0) GARCH <- TRUE
    if( (p == 0) && (q == 0) ) {ARMAonly = TRUE}
    optim.finished <- FALSE
    if (AR == TRUE)
        m <- 1
    if( MA == TRUE) 
        n <- 1
    if (GARCH == TRUE)
        q <- 1
    mn <- max(m,n); pq <- max(p,q)
    garchLLH = function(parm){
        # model parameters
        if(sum(is.nan(parm)) != 0) {return(1e99)}
        mu <- parm[1];
        a <- parm[(1+1):(2+m-1)]; b <- parm[(1+m+1):(2+m+n-1)]
        omega <- parm[1+m+n+1]; alpha <- parm[(2+m+n+1):(3+m+n+p-1)]
        gm <- parm[(2+m+n+p+1):(3+m+n+p+p-1)]
        beta <- parm[(2+m+n+2*p+1):(3+m+n+2*p+q-1)]
        delta <- parm[2+m+n+2*p+q+1]; 
        skew <- parm[3+m+n+2*p+q+1]; shape <- parm[4+m+n+2*p+q+1];
        if( !APARCH ) 
        { 
            gm = rep(0,p);
            delta = 2; 
            if( cond.dist == "stable" ) delta = 1
        }
        # Setting parameters to accept tapper off MA, AR or GARCH coefficients
        if( AR == TRUE) 
        	a <- 0
        if( MA == TRUE) 
        	b <- 0
        if( GARCH == TRUE) 
        	beta <- 0
        if (intercept == FALSE)
            mu <- 0
	  # Avoid being out of parameter space
        parset <- c(omega,alpha,if(!GARCH) beta,delta)
        cond.general <- any(parset < GStol) 
        cond.normal <- FALSE
        cond.student <- FALSE
        cond.gev <- FALSE
        cond.stable <- FALSE
        if( cond.dist == "norm")
           cond.normal <- ( sum(alpha) + sum(beta) > 1 - GStol )
        if( cond.dist == "stable")
        {
            if( shape-delta < GSstable.tol || abs(shape) < GSstable.tol ||
                !(abs(skew) < 1) || !((shape - 2) < 0) )
            {
                return(1e99)
            }
            tau <- skew*tan(shape*pi/2)
            kdelta <- pi/2
            if(abs(delta-1) > GStol) kdelta <- gamma(1 - delta)*cos(pi*delta/2)
        	lamb <- kdelta^(-1)*gamma(1 - delta/shape)*(1 + tau^2)^(delta/2/shape)*
					cos(delta/shape*atan(tau))
            cond.stable <- FALSE
        }
        if (cond.general || cond.student || cond.gev || cond.stable)
        {
            return(1e99)
        }

        # tapper off the mean equation 2
        e.init <- rep(0,mn)
        e.parc <- c(e.init,filter(data, filter = c(1, -a), sides = 1)[(mn+1):N])
        e.res <- c( e.init, filter(e.parc[-(1:mn)], filter = -b,
        method = "recursive", init = e.init[1:n]))     

        # find i.i.d sequence z
        z <- e.res - mu
        if(ARMAonly)
            hh <- rep(omega,N)
        else
        {
            h <- rep(0.1, pq)
            edeltat = 0
            for( i in 1:p)
            {
                edelta <- alpha[i]*(abs(z)-gm[i]*z)^delta
                edeltat = edeltat +  edelta[(p-(i-1)):(N-i)]
            }
            edeltat = c(h[1:p],edeltat)
            c <- omega/(1-sum(beta))
            h <- c( h[1:pq], c + filter(edeltat[-(1:pq)], filter = beta,
            method = "recursive", init = h[q:1]-c))
            hh <- abs(h)^(1/delta)
        }

        # get output Residuals ?
        if (optim.finished & get.res)
        {
           out$ARMA.res <<- z
           out$GARCH.sig <<- hh
	  }
 
        # Return llh function        
        llh.dens <- GSgarch.Dist(z = z, hh = hh, shape = shape, 
        skew = skew, cond.dist = cond.dist)
        llh <- llh.dens
        if (is.nan(llh) || is.infinite(llh) ||
            is.na(llh)) 
        {
            llh <- 1e99
        }
        llh
    }
    # Performing optimization
    start <- GSgarch.GetStart(data = data,m = m,n = n,p = p,q = q,AR = AR,
             MA = MA, cond.dist = cond.dist)
    rest <- function(parm)
    {
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
    if (algorithm == "sqp")
        fit1 <- solnp(pars = start[1,], fun = garchLLH, 
                LB = start[2,], UB = start[3,], control = control)
    if (algorithm == "sqp.rest")
        fit1 <- solnp(pars = start[1,], fun = garchLLH, ineqfun = rest, ineqLB = 0,
                ineqUB = 1, LB = start[2,], UB = start[3,], control = control)
    if (algorithm == "nlminb")
    {
        fit1 <- nlminb(start[1,], objective = garchLLH,
	          lower = start[2,], upper = start[3,], 
	          control = control)
        out$llh <- fit1$objective
        out$par <- fit1$par
        out$hessian <- optim(par = fit1$par, fn = garchLLH, 
	                 method = "Nelder-Mead", hessian = TRUE)$hessian 
    }
    if (any(c("sqp", "sqp.rest") == algorithm))
    {
        out$llh <- fit1$values[length(fit1$values)]
        out$par <- fit1$pars
        out$hessian <- fit1$hessian
    } 

    # Organizing the output of the program
    optim.finished = TRUE
    if(get.res)
       garchLLH(out$par)
    outindex <-   c(if(intercept) 1, 
               if(AR == FALSE) (1+1):(2+m-1),
               if(MA == FALSE) (1+m+1):(2+m+n-1),
               (1+m+n+1),
               if(!ARMAonly) (2+m+n+1):(3+m+n+p-1),
               if(APARCH) (2+m+n+p+1):(3+m+n+p+p-1),
               if(!GARCH) (2+m+n+2*p+1):(3+m+n+2*p+q-1),
               if(APARCH) (2+m+n+2*p+q+1),
               if(any(c("sstd","stable")  == cond.dist)) (3+m+n+2*p+q+1),
               if(any(c("std","gev","stable","sstd")  == cond.dist)) 
               (4+m+n+2*p+q+1))
    outnames <- c(if(intercept) "mu", if( AR == FALSE) paste("ar", 1:m, sep = ""),
               if(MA == FALSE) paste("ma", 1:n, sep = ""),
               "omega",
               if(!ARMAonly) paste("alpha", 1:p, sep = ""),
               if(APARCH) paste("gamma", 1:p, sep = ""),
               if(!GARCH) paste("beta", 1:q, sep = ""),
               if(APARCH) "delta",
               if(any(c("sstd","stable")  == cond.dist)) "skew",
               if(any(c("std","gev","stable","sstd")  == cond.dist)) 
               "shape")
    out$par <- out$par[outindex]
    names(out$par) <- outnames
    out$hessian <- out$hessian[outindex,outindex]
    nParam <- length(out$par)
    out$aic  = 2*out$llh + 2*nParam 
    out$aicc = 2*out$llh + 2*nParam*N/(N - nParam - 1)
    out$bic =  2*out$llh + nParam*log(N)

    # Print Summary
    if ( printRes ) 
    {
		solveHessianFailed = FALSE
		out$se.coef <- 0
        out$se.coef <- try(sqrt(diag(solve(out$hessian))), silent = TRUE)
		if(!is.numeric(out$se.coef))
		{
			solveHessianFailed = TRUE
			print("Error solving Hessian Matrix. Variable 'se.coef' have the reported error.")
			out$matcoef <- cbind(out$par, rep(NA,length(out$par)), rep(NA,length(out$par)), rep(NA,length(out$par)))
			dimnames(out$matcoef) = dimnames(out$matcoef) = list(names(out$par), c(" Estimate",
			" Std. Error", " t value", "Pr(>|t|)"))
		}
		else
		{
			out$tval <- try(out$par/out$se.coef, silent = TRUE)
			out$matcoef = cbind(out$par, if(is.numeric(out$se.coef)) out$se.coef, if(is.numeric(out$tval)) out$tval, 
                      if(is.numeric(out$tval)) 2*(1-pnorm(abs(out$tval))))
			dimnames(out$matcoef) = list(names(out$tval), c(" Estimate",
			" Std. Error", " t value", "Pr(>|t|)"))
		}
		cat("\nFinal Estimate of the Negative LLH:\n")
		cat("LLH:",out$llh)
		if(out$llh == 1e99)
			print("Algorithm did not achieved convergence.")
		cat("\nCoefficient(s):\n")
		printCoefmat(round(out$matcoef,digits=6), digits = 6, signif.stars = TRUE)
    }
    return(out)
}
############
# Simulation 
############
GSgarch.Sim <- function(N = 1000,mu = 0.1,a = c(0.5,0.3),
b = c(-0.4,0.3,-0.1),omega = 0.05, alpha = c(0.1),gm = c(0),
beta = c(0.1,0.05,0.03), delta = 2, skew = 0, shape = 3, cond.dist = "norm")
{
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
    # we will not verify stationarity conditions for simulation to 
    # allow the user simulate also non stationary series
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
##############
# AIC fit
##############
GSgarch.GetOrder <- function(m,n,p,q)
{
	# error treatment on input parameters
    if(m%%1 != 0 || n%%1 != 0 || p%%1 != 0 || q%%1 != 0 || 
       any (c(m,n,p,q) < 0) || (p == 0 && q != 0) ||
       any (c(m,n,p,q) > 10) ) 
       stop ("Invalid ARMA-GARCH order. We allow pure GARCH or APARCH. AR/MA/ARMA-GARCH/APARCH models.
	          The order of the parameters could be set up to 10.")
    arma.garch.order <- c()
    for(i1 in 0:m)
    {
        for(i2 in 0:n)
        {
            for(i3 in 1:p)
            {
                for(i4 in 0:q)
                { 
                   ord <- c(i1,i2,i3,i4)
                   arma.garch.order <- rbind(arma.garch.order,c(i1,i2,i3,i4))
                }
            }
        }
    }
    return(arma.garch.order)
}
GSgarch.FitAIC <- function(data,mMAX=1,nMAX=1,pMAX=1,qMAX=1, cond.dist = "norm", 
                  algorithm = "sqp",APARCH = FALSE, intercept = TRUE,control = NULL)
{
	# error treatment on input parameters
    cond.dist.list <- c("norm", "std", "sstd", "gev", "stable")
    if( !any(cond.dist.list == cond.dist) )   
        stop ("Invalid Conditional Distribution. Choose: norm,std,sstd,gev or stable")
	if( !is.numeric(data) || !is.vector(data))
		stop("data set must be a numerical one dimensional vector")
		
	# begin of function
    T <- length(data)
    aic.min <- 1e10
    aic <- 1e10
    fit.min <- list()
    fit <- list()
    aic.list <- GSgarch.GetOrder(mMAX,nMAX,pMAX,qMAX)
    aic.list.size <- length(aic.list[,1])
    for( i in 1:aic.list.size)
    {
        fit = GSgarch.Fit(data, m = aic.list[i,1],n = aic.list[i,2],
              p = aic.list[i,3],q = aic.list[i,4], cond.dist = cond.dist, 
              APARCH = APARCH, algorithm = algorithm, 
              intercept = intercept, control = control)
        nParam <- length(fit$par)
        aic  = 2*fit$llh + 2*nParam 
        aicc = 2*fit$llh + 2*nParam*T/(T - nParam - 1)
        bic =  2*fit$llh + nParam*log(T)
        cat(aic.list[i,],"-llh:",fit$llh,"AIC:",aic," AICC:",aicc," BIC:",bic,"\n")
        if (aic < aic.min)
        {
            fit.min <- fit
            fit.min$order <- aic.list[i,]
            aic.min <- aic
        }
    }
    return(fit.min)
}