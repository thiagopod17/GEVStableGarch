
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
#  GSgarch.Fit             Fits ARMA-GARCH or ARMA-APARCH model  
################################################################################

GSgarch.Fit <-
function(
    formula = ~ garch(1,1), 
    data,  
    cond.dist = c("norm", "std", "sstd", "gev", "stable"),
    include.mean = TRUE, 
    algorithm = c("sqp","sqp.restriction","nlminb"),
    printRes = TRUE,
    control = NULL,
    DEBUG = FALSE)
{  
    # Description:
    #     This functions reads the univariate time series and fits
    #     a ARMA-GARCH/APARCH model with conditional GEV and stable
    #     distribution.
    #     TIME SERIES MODEL:
    #         According to Wurtz et al. (2006) 
    #         xt = mu + a(B)xt + b(B)et,
    #         et = zt * sigmat
    #         zt ~ Dv(0, 1)
    #         sigmat^2 = omega + sum(alphai(e(t-i))^2) + sum(betaj*sigmat-j^2)
    #     REMARKS:
    #     et for ARMA(m,n), n > 1 will be initiated as 0, i.e, e[-n+1:0] = 0.
    #     zt will be initiated as 0.
    #     ht will be initiated as 0.1 (see eq. (22) of Wurtz, 2006).   
    #     VARIABLE NOTATION USED INSIDE THIS FUNCTION:
    #         N: sample size
    #         m,n: ARMA order
    #         p,q, pq = max(p,q): GARCH order. They were called u,v, uv in Wurtz et al.(2006)
    #         mu, a, b: ARMA parameters
    #         omega, alpha, gamma, beta: APARCH vector parameters
    #         garchLLH: Log Likelihood for the ARMA(m,n)-GARCH(p,q) model
    #         parm: ARMA-GARCH parameters concatenated as c(mu,a,b,omega,alpha,gamma,beta)
    #         h: conditional variance. It was referenced as sigmat in model equations
    #         x: data set. It is the time series that will be modelled as ARMA-GARCH
    #         xi: GEV shape parameter. xi > -0.5 (llh) and < 0.5 (finiteness of variance, see Zhao et al. (2011))
    #         AR, MA: if AR (or MA) equals TRUE, then the model has AR (or MA) order equals to zero.
    #         param or pm (stabledist): Parametrization of stable distributions as chosen to be 2.  
    #         Even for stable distribution llh, there's a problem for finding the estimators for
    #         alpha near to 2 and beta > 0. The bad performance on the ARMA-GARCH model was
    #         similar to the single llh estimation of i.i.d samples of stable distribution. 
    #         In our simulations, the param = 2 was chosen to be the ones that performs
    #         better in these situations.
    #     IMPORTANT DETAILS FOR USING THIS FUNCTION:
    #         Some care must be taken when using this function to estimate the parameters
    #         of some models. 
    #         For pure GARCH(p,q) with p,q >= 1 we need to make sure that the 'gamma' variable
    #         is a vector of length 'p' and with all entries equal to zero.
    #         For pure ARCH(p) = GARCH(p,0) we set the variable GARCH equal to TRUE to indicate
    #         that the model has order q = 0. Then, we make q = 1 to estimate the parameters 
    #         of a GARCH(p,1) with beta = 0 and gamma = 0. 
         
    # Arguments:
    #   formula - ARMA(m,n) + GARCH/APARCH(p,q) mean and variance specification 
    #   data - vector of data
    #   m, n, p, q - model order as in ARMA(m,n)-GARCH/APARCH(p,q)
    #   include.mean - a logical, should the mean value be estimated ? 
    #   algorithm - 
    #   cond.dist - name of the conditional distribution, one of
    #       gev, stable, norm, std, sstd   
    
    # Return:
    #   Asdf - The asdf     
      
    # FUNCTION:  

    # Error Treatment on input parameters
    cond.dist = match.arg(cond.dist)
    algorithm = match.arg(algorithm)
    if (!is.numeric(data) || any(is.na(data)) || any(is.null(data)) || any(!is.finite(data)))
        stop("Invalid 'data' input. It may be contain NA, NULL or Inf.")     
      
    # Call:
    CALL = match.call()  
  
    # Configuring Tolerance
    GSstable.tol = 1e-2 # upper and lower bounds tolerance. Should be greater than tol
    GStol = 1e-8 # General tolerance for arma-garch parameters. 
    # In the beggining it was set to 1e-5

    # Getting order model from object formula
    formula <- .getFormula(formula)
    m <- formula$formula.order[1]
    n <- formula$formula.order[2]    
    p <- formula$formula.order[3]
    q <- formula$formula.order[4]
    APARCH <- formula$isAPARCH
        
    # Configuring model order
    ARMAonly <- FALSE
    AR <- FALSE 
    MA <- FALSE 
    GARCH <- FALSE
    if( m == 0) AR <- TRUE
    if( n == 0) MA <- TRUE
    if( q == 0) GARCH <- TRUE
    optim.finished <- FALSE
    if (AR == TRUE)
        m <- 1
    if( MA == TRUE) 
        n <- 1
    if( (p == 0) && (q == 0))
        ARMAonly = TRUE
    if (GARCH == TRUE && !ARMAonly)
        q <- 1
    
    # Initial configurations
    data <- data; 
    N <- length(data)
    out <- NULL # output of the GSgarch.Fit function
    out$order <- c(m,n,p,q,include.mean,APARCH)
    TMPvector <- c(if(m != 0 || n != 0) c("arma(",m,",",n,")-"),
                   if(APARCH==FALSE)c("garch(",p,",",q,")"),
                   if(APARCH==TRUE)c("aparch(",p,",",q,")"))
    TMPorder <- paste(TMPvector, collapse="")
    TMPvectorintercept <- c("include.mean:",if(include.mean==TRUE)"TRUE", 
                            if(include.mean==FALSE)"FALSE")
    TMPintercept <- paste(TMPvectorintercept, collapse="")
    out$model <- paste(TMPorder,"##",TMPintercept, collapse="")
    out$cond.dist <- cond.dist
    out$data <- data
    optim.finished <- FALSE
    

    
    # BEGIN: Log Likelihood for pure ARMA process
    #############################################
    armaLLH <- function(parm)
    {      
        if(DEBUG)
            print(parm)
        # check if some parameters are NAN
        if(sum(is.nan(parm)) != 0) {return(1e99)}
        
        # Getting parameters from parm vector 
        mu <- parm[1];
        a <- parm[(1+1):(2+m-1)]
        b <- parm[(1+m+1):(2+m+n-1)]
        skew <- parm[1+m+n+1]
        shape <- parm[2+m+n+1]
        
        # Setting parameters to accept tapper off MA, AR or GARCH coefficients
        if( AR == TRUE) 
          a <- 0
        if( MA == TRUE) 
          b <- 0
        if (include.mean == FALSE)
          mu <- 0
         
        if( cond.dist == "stable")
        {
          if( shape-delta < GSstable.tol || abs(shape) < GSstable.tol ||
                !(abs(skew) < 1) || !((shape - 2) < 0) )
          {
            return(1e99)
          }
        }
        
        # Filters the Time series to obtain the i.i.d. sequence of 
        # 'innovations' to evaluate the Log-likelihood function
        z <- filter.Arma(data = data, m = m, n = n, mu = mu, a = a, b = b)
        
        if(DEBUG)
        {
            print(c("length(z)",length(z)))
            print(c("N",N))
        }
        
        
        # get output Residuals
        if (optim.finished)
        {
          out$ARMA.res <<- z
          out$GARCH.sig <<- 1
        }
        
        # Return llh function        
        llh.dens <- GSgarch.Dist(z = z, hh = 1, shape = shape, 
                                 skew = skew, cond.dist = cond.dist)
        llh <- llh.dens
        if (is.nan(llh) || is.infinite(llh) || is.na(llh)) 
        {
          llh <- 1e99
        }
        llh
    }
    # END: Log Likelihood for pure ARMA process
    #############################################


    # BEGIN: garch Likelihood ARMA-APARCH or pure GARCH process
    ###########################################################
    garchLLH = function(parm){
        
        # check if some parameters are NAN
        if(sum(is.nan(parm)) != 0) {return(1e99)}
        
        # Getting parameters from parm vector 
        mu <- parm[1];
        a <- parm[(1+1):(2+m-1)]
        b <- parm[(1+m+1):(2+m+n-1)]
        omega <- parm[1+m+n+1]
        alpha <- parm[(2+m+n+1):(3+m+n+p-1)]
        gm <- parm[(2+m+n+p+1):(3+m+n+p+p-1)]
        beta <- parm[(2+m+n+2*p+1):(3+m+n+2*p+q-1)]
        delta <- parm[2+m+n+2*p+q+1] 
        skew <- parm[3+m+n+2*p+q+1]
        shape <- parm[4+m+n+2*p+q+1]
        
        
        # Configuring delta and gamma for Garch estimation
        if( !APARCH) 
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
        if (include.mean == FALSE)
            mu <- 0
        
        # Avoid being out of parameter space
        cond.general <- FALSE
        cond.normal <- FALSE
        cond.student <- FALSE
        cond.gev <- FALSE
        cond.stable <- FALSE
        parset <- c(omega,alpha,if(!GARCH) beta,delta)
        cond.general <- any(parset < GStol) 
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
            if(abs(delta-1) > GStol) 
                kdelta <- gamma(1 - delta)*cos(pi*delta/2)
                lamb <- kdelta^(-1)*gamma(1 - delta/shape)*(1 + tau^2)^(delta/2/shape)*
                    cos(delta/shape*atan(tau))
            cond.stable <- FALSE
        }

        # out of parameter space
        if (cond.general || cond.student || cond.gev || cond.stable)
        {
            return(1e99)
        }
        
        # Filters the Time series to obtain the i.i.d. sequence of 
        # 'innovations' to evaluate the Log-likelihood function
        if(AR == TRUE && MA  == TRUE)
        {
            filteredSeries <- filter.Aparch(data = data,p = p,q = q, 
              mu = mu, omega = omega, alpha = alpha, beta = beta, gamma = gm, delta = delta)
            z <- filteredSeries[,1]
            hh <- filteredSeries[,2]          
        }
        if(AR == FALSE || MA == FALSE)
        {
            filteredArma <- filter.Arma(data = data, m = m, n = n, mu = mu, a = a, b = b)
            filteredSeries <- filter.Aparch(data = filteredArma,p = p,q = q, 
                              mu = 0, omega = omega, alpha = alpha, beta = beta, gamma = gm, delta = delta)
            z <- filteredSeries[,1]
            hh <- filteredSeries[,2]          
          
        }
          # CURRENT FILTERING THAT WORKS REALLY GOOD FOR PURE APARCH PROCESS
        
#         # if only garch(p,0) or aparch(p,0)
#         if(GARCH == TRUE)
#           beta = 0
#         Mean.z <- mean(abs(z)^delta)
#         for( i in 1:p)
#         {
#           edelta <- alpha[i]*(c(rep(Mean.z,p),((abs(z)-gm[i]*z)^delta)[1:(N-1)]))
#           edeltat = edeltat +  edelta[(p-(i-1)):(p+N-i)]
#         }
#         edeltat = omega + edeltat
#         
#         h <- filter(edeltat, filter = beta,
#                     method = "recursive", init = rep(Mean.z,q))
#         hh <- abs(h)^(1/delta)



        # get output Residuals
        if (optim.finished)
        {
            out$ARMA.res <<- z
            out$GARCH.sig <<- hh
        }
        
        # Return llh function        
        llh.dens <- GSgarch.Dist(z = z, hh = hh, shape = shape, 
                                 skew = skew, cond.dist = cond.dist)
        llh <- llh.dens
        if (is.nan(llh) || is.infinite(llh) || is.na(llh)) 
        {
            llh <- 1e99
        }
        llh
    }
    # END: garch Likelihood ARMA-APARCH or pure GARCH process
    ###########################################################


    # Getting start, lower and upper bound for parameters to perform optimization
    start <- GSgarch.GetStart(data = data,m = m,n = n,p = p,q = q,AR = AR,
                              MA = MA, ARMAonly = ARMAonly, cond.dist = cond.dist)
    if(DEBUG)
    {
        print("start")
        print(start)
    }
    
    # Function that evaluate stationarity conditions to guide parameter estimation.
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
    if (algorithm == "sqp.restriction")
        fit1 <- solnp(pars = start[1,], fun = garchLLH, ineqfun = rest, ineqLB = 0,
                    ineqUB = 1, LB = start[2,], UB = start[3,], control = control)
    if (algorithm == "nlminb")
    {      
        if(ARMAonly)
        {            
          fit1 <- nlminb(start[1,], objective = armaLLH,
                       lower = start[2,], upper = start[3,], 
                       control = control)
          out$llh <- fit1$objective
          out$par <- fit1$par
          out$hessian <- optim(par = fit1$par, fn = armaLLH, 
                               method = "Nelder-Mead", hessian = TRUE)$hessian
        }
        else{
            fit1 <- nlminb(start[1,], objective = garchLLH,
                       lower = start[2,], upper = start[3,], 
                       control = control)
            out$llh <- fit1$objective
            out$par <- fit1$par
            out$hessian <- optim(par = fit1$par, fn = garchLLH, 
                                 method = "Nelder-Mead", hessian = TRUE)$hessian 
        }
    }
    if (any(c("sqp", "sqp.restriction") == algorithm))
    {
        out$llh <- fit1$values[length(fit1$values)]
        out$par <- fit1$pars
        out$hessian <- fit1$hessian
    } 
    
    if(DEBUG)
      print(fit1)

    # Organizing the output of the program
    optim.finished = TRUE
    
    # Call garchLLH function to update the values of the ARMA residuals 
    # and the GARCH/APARCH volatility.
    if(ARMAonly)        
        armaLLH(out$par)
    else 
        garchLLH(out$par)
    
    # Creating index to create a vector with the estimated parameters.
    if(!ARMAonly)
    {
        outindex <- c(if(include.mean) 1, 
                    if(AR == FALSE) (1+1):(2+m-1),
                    if(MA == FALSE) (1+m+1):(2+m+n-1),
                    if(!ARMAonly) (1+m+n+1),
                    if(!ARMAonly) (2+m+n+1):(3+m+n+p-1),
                    if(APARCH) (2+m+n+p+1):(3+m+n+p+p-1),
                    if(!GARCH) (2+m+n+2*p+1):(3+m+n+2*p+q-1),
                    if(APARCH) (2+m+n+2*p+q+1),
                    if(any(c("sstd","stable")  == cond.dist)) (3+m+n+2*p+q+1),
                    if(any(c("std","gev","stable","sstd")  == cond.dist)) 
                      (4+m+n+2*p+q+1))
    } else {
        outindex <- c(if(include.mean) 1, 
                    if(AR == FALSE) (1+1):(2+m-1),
                    if(MA == FALSE) (1+m+1):(2+m+n-1),
                    if(!ARMAonly) (1+m+n+1),
                    if(!ARMAonly) (2+m+n+1):(3+m+n+p-1),
                    if(APARCH) (2+m+n+p+1):(3+m+n+p+p-1),
                    if(!GARCH) (2+m+n+2*p+1):(3+m+n+2*p+q-1),
                    if(APARCH) (2+m+n+2*p+q+1),
                    if(any(c("sstd","stable")  == cond.dist)) (1+m+n+2*p+q+1),
                    if(any(c("std","gev","stable","sstd")  == cond.dist)) 
                      (2+m+n+2*p+q+1))  
    }
    
    outnames <- c(if(include.mean) "mu", 
                  if( AR == FALSE) paste("ar", 1:m, sep = ""),
                  if(MA == FALSE) paste("ma", 1:n, sep = ""),
                  if(!ARMAonly) "omega",
                  if(!ARMAonly) paste("alpha", 1:p, sep = ""),
                  if(APARCH) paste("gamma", 1:p, sep = ""),
                  if(!GARCH) paste("beta", 1:q, sep = ""),
                  if(APARCH) "delta",
                  if(any(c("sstd","stable")  == cond.dist)) "skew",
                  if(any(c("std","gev","stable","sstd")  == cond.dist)) 
                    "shape")
    
    if(DEBUG)
    {
        print(c("out",out))
        print(c("outindex",outindex))
    }
      
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
            dimnames(out$matcoef) = dimnames(out$matcoef) = 
            list(names(out$par), c(" Estimate"," Std. Error", " t value", "Pr(>|t|)"))
        }
        else
        {
            out$tval <- try(out$par/out$se.coef, silent = TRUE)
            out$matcoef = cbind(out$par, if(is.numeric(out$se.coef)) out$se.coef, if(is.numeric(out$tval)) out$tval, 
                                if(is.numeric(out$tval)) 2*(1-pnorm(abs(out$tval))))
            dimnames(out$matcoef) = list(names(out$tval), 
                                    c(" Estimate"," Std. Error", " t value", "Pr(>|t|)"))
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



################################################################################

