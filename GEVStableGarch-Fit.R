
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

garchFit
GSgarch.Fit <-
function(
    formula = ~ garch(1,1), 
    data,  
    cond.dist = c("norm", "std", "sstd", "gev", "stable"),
    include.mean = TRUE, 
    algorithm = c("sqp","sqp.restriction","nlminb"),
    printRes = TRUE,
    control = NULL)
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
         
    # Arguments:
    #   formula - ARMA(m,n) + GARCH/APARCH(p,q) mean and variance specification 
    #   data - vector of data
    #   m, n, p, q - model order as in ARMA(m,n)-GARCH/APARCH(p,q)
    #   include.mean - a logical, should the mean value be estimated ? 
    #   algorithm - 
    #   cond.dist - name of the conditional distribution, one of
    #       gev, stable, norm, std, sstd   
    #   GStol.b - upper and lower bounds tolerance. Should be greater than tol
    #   GStol - General tolerance for arma-garch parameters. In the beggining it was set to 1e-5    
    
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
    GSstable.tol = 1e-2 
    GStol = 1e-8

    # Getting order model from object formula
    formula <- .getFormula(formula)
    m <- formula$formula.order[1]
    n <- formula$formula.order[2]    
    p <- formula$formula.order[3]
    q <- formula$formula.order[4]
    APARCH <- formula$isAPARCH
    
#     # Checking if model order was specified correctly
#     if(m%%1 != 0 || n%%1 != 0 || p%%1 != 0 || q%%1 != 0 || 
#         any (c(m,n,p,q) < 0) || (p == 0 && q != 0) || (p == 0 && APARCH) ||
#         any (c(m,n,p,q) > 10) ) 
#         stop ("Invalid ARMA-GARCH order. We allow pure GARCH or APARCH. AR/MA/ARMA-GARCH/APARCH models.
#             The order of the parameters could be set up to 10.")

    # Initial configurations
    data <- data; N <- length(data)
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
    ARMAonly <- FALSE
    if( (p == 0) && (q == 0))
    {
        ARMAonly = TRUE
        p = 1
    }
    optim.finished <- FALSE
    
    # Configuring model order
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
        if (include.mean == FALSE)
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

        # out of parameter space
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
    
    # Performing optimization
    start <- GSgarch.GetStart(data = data,m = m,n = n,p = p,q = q,AR = AR,
                              MA = MA, cond.dist = cond.dist)
    
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
        fit1 <- nlminb(start[1,], objective = garchLLH,
                       lower = start[2,], upper = start[3,], 
                       control = control)
        out$llh <- fit1$objective
        out$par <- fit1$par
        out$hessian <- optim(par = fit1$par, fn = garchLLH, 
                             method = "Nelder-Mead", hessian = TRUE)$hessian 
    }
    if (any(c("sqp", "sqp.restriction") == algorithm))
    {
        out$llh <- fit1$values[length(fit1$values)]
        out$par <- fit1$pars
        out$hessian <- fit1$hessian
    } 
    
    # Organizing the output of the program
    optim.finished = TRUE
    
    # Call garchLLH function to update the values of the ARMA residuals 
    # and the GARCH/APARCH volatility.
    garchLLH(out$par)
    
    # Creating index to create a vector with the estimated parameters.
    outindex <-   c(if(include.mean) 1, 
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
    
    outnames <- c(if(include.mean) "mu", if( AR == FALSE) paste("ar", 1:m, sep = ""),
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







################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .getFormula             Gets the formula.mean and formula.variance from
#                          the object formula.
################################################################################


.getFormula <-
function(
    formula)
{
    # Description:
    #   This functions reads formula object and converts it into a list 
    #   containing the separeted mean and variance arguments.
    #   Examples from output:
    #     ~arma(1,1)+garch(1,1): formula.mean = ~arma(1,1); formula.variance = ~garch(1,1)
    #     ~aparch(1,1): formula.mean = ~arma(0,0); formula.variance = ~aparch(1,1)
    #     ~arch(1): formula.mean = ~arma(0,0); formula.variance = ~arch(1)
    #     ~arma(1,1): formula.mean = ~arma(1,1); formula.variance = ~garch(0,0)
    #     ~ar(1): formula.mean = ~ar(1); formula.variance = ~garch(0,0) 
    #     ~ma(1): formula.mean = ~ma(1); formula.variance = ~garch(0,0) 
  
    # Arguments:
    #   formula - ARMA(m,n) + GARCH/APARCH(p,q) mean and variance specification 
    
    # Return:
    #   A list containing two elements, formula.mean and formula.variance     
    
    # FUNCTION: 
  
    # Initial variable declaration
    allLabels = attr(terms(formula), "term.labels")
    formulas.mean.allowed = c("arma")
    formulas.variance.allowed = c("garch","aparch")
    formulaOK <- TRUE
    checkFormulaMean <- ""
    checkFormulaVariance <- ""
    isAPARCH = FALSE
    
    # Error treatment of input parameters 
    if( (length(allLabels) != 1 ) && (length(allLabels) != 2 ) )
        formulaOK <- FALSE
    
    # Formula of type: ~ formula1 + formula2
    else if (length(allLabels) == 2)
    {
        formula.mean = as.formula(paste("~", allLabels[1]))
        formula.var = as.formula(paste("~", allLabels[2]))
        checkFormulaMean = rev(all.names(formula.mean))[1]
        checkFormulaVariance = rev(all.names(formula.var))[1]
        if( !any(formulas.mean.allowed == checkFormulaMean) || 
            !any(formulas.variance.allowed == checkFormulaVariance))   
            formulaOK <- FALSE
    }
    # Formula of type: ~formula1
    else if (length(allLabels) == 1) 
    {
      
        # pure 'garch' or 'aparch'
        if(grepl("arch", attr(terms(formula), "term.labels"))) 
        {
            formula.mean = as.formula("~ arma(0, 0)")
            formula.var = as.formula(paste("~", allLabels[1]))
            checkFormulaVariance = rev(all.names(formula.var))[1]
            if(!any(formulas.variance.allowed == checkFormulaVariance))   
                formulaOK <- FALSE
        }
        else # pure 'ar', 'ma' or 'arma' model.
        {
            formula.mean = as.formula(paste("~", allLabels[1]))  
            formula.var = as.formula("~ garch(0, 0)") 
            checkFormulaMean = rev(all.names(formula.mean))[1]
            if(!any(formulas.mean.allowed == checkFormulaMean))   
                formulaOK <- FALSE
        }    
    }

    # Check if we are fitting "aparch" model 
    if(checkFormulaVariance == "aparch")
        isAPARCH = TRUE
    
    # Get model order and check if they were specified correctly
    if(formulaOK == TRUE)
    {
        model.order.mean = 
            as.numeric(strsplit(strsplit(strsplit(as.character(formula.mean), 
            "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
        model.order.var = 
            as.numeric(strsplit(strsplit(strsplit(as.character(formula.var), 
            "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
        if( (length(model.order.mean) != 2) || (length(model.order.mean) != 2))
          formulaOK <- FALSE 
    }
    
    # Check if model order was specified correctly.
    if(formulaOK == TRUE)
    {
        m = model.order.mean[1]
        n = model.order.mean[2]
        p = model.order.var[1]
        q = model.order.var[2]        
        if(m%%1 != 0 || n%%1 != 0 || p%%1 != 0 || q%%1 != 0 || 
            any (c(m,n,p,q) < 0) || (p == 0 && q != 0))
            formulaOK <- FALSE            
    }
    
    
    # Stop if formula was not specified correctly
    if(formulaOK == FALSE)
        stop ("Invalid Formula especification. 
            Formula mean must be 'arma' and 
            Formula Variance must be one of: garch or aparch
            For example:
                ARMA(1,1)-GARCH(1,1):  ~arma(1,1)+garch(1,1),
                AR(1)-GARCH(1,1):      ~arma(1,0)+garch(1,1),
                MA(1)-APARCH(1,0):     ~arma(0,1)+aparch(1,0),
                ARMA(1,1):             ~arma(1,1),
                ARCH(2):               ~garch(1,0),
            For more details just type: ?GSgarch.Fit")
    
    # Return
    list(formula.mean = formula.mean,formula.var = formula.var, 
         formula.order = c(m,n,p,q), isAPARCH = isAPARCH)
}






