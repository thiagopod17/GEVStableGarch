#' Generalized Asymmetric t (GAT) Distribution
#'
#' Functions to compute the density, distribution function, quantile function,
#' and to generate random variates for the Generalized Asymmetric t (GAT)
#' distribution.
#'
#' This distribution corresponds to the \code{t3} distribution described in
#' Paolella (1997) and Mittnik and Paolella (2000). The GAT family includes
#' the Student's t, Laplace, Cauchy, and normal distributions as special cases.
#' In particular, as \eqn{\nu \to \infty} the distribution approaches normality.
#'
#' @name GAT
#' @rdname GAT
#' @aliases gat dgat pgat qgat rgat
#'
#' @param x Numeric vector of values for calculating density. 
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations for random generation.
#' @param mean Location parameter.
#' @param sd Scale parameter (must be > 0).
#' @param nu Tail (shape) parameter (must be > 0).
#' @param d Second shape parameter (must be > 0).
#' @param xi Asymmetry parameter (must be > 0; \eqn{\xi = 1} gives symmetry).
#' @param log Logical; if \code{TRUE}, densities are returned on the log scale.
#'
#' @return
#' \item{dgat}{density values}
#' \item{pgat}{distribution function values}
#' \item{qgat}{quantile function values}
#' \item{rgat}{random variates}
#'
#' @references
#' Mittnik, S., Paolella, M. S. (2000).
#' Prediction of Financial Downside-Risk with Heavy-Tailed Conditional
#' Distributions.
#'
#' Paolella, M. (1997).
#' Tail Estimation and Conditional Modeling of Heteroskedastic Time-Series.
#' PhD Thesis, Institute of Statistics and Econometrics,
#' Christian Albrechts University of Kiel.
#' 
#' Bertocchi, M., Giacometti, R., Ortobelli, S., & Rachev, S. T. (2005).
#' The impact of different distributional hypothesis on returns in asset allocation.
#' \emph{Finance Letters}, 3(1), 17-27.
#'
#' @author Thiago do Rego Sousa
#'
#' @examples
#' par(mfrow = c(2, 2))
#' set.seed(1000)
#' r <- rgat(n = 1000)
#' plot(r, type = "l", main = "GAt Random Values")
#'
#' hist(r, probability = TRUE, border = "white")
#' x <- seq(min(r), max(r), length = 201)
#' lines(x, dgat(x), lwd = 2)
#'
#' plot(sort(r), (1:1000)/1000, main = "Probability", ylab = "Probability")
#' lines(x, pgat(x), lwd = 2)
#'
#' round(qgat(pgat(q = seq(-10, 10, by = 0.5))), 6)
#'
#' @export
dgat <- 
  function(x, mean = 0, sd = 1, nu = 2, d = 3, xi = 1, log = FALSE)
  {   

  
    
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


#' @export
pgat <- 
  function(q, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
  {   
    
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
        return ( replace(pw / ( nu + pw ),which (pw == Inf, arr.ind = TRUE),1) )
    } 

    # Compute auxiliary variables:
    z = (q - mean ) / sd
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


#' @export
qgat <- 
  function(p, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)  
  {   
    
    # Define auxiliary functions
    Lp <- function (p = p, nu = nu, d = d, xi = xi)
    {
        qbeta( ( 1 + xi^2 ) * p, nu, 1/d)
    }
    Up <- function (p = p, nu = nu, d = d, xi = xi)
    {
        qbeta( ( p - 1 / ( 1 + xi^2 ) ) * ( 1 + xi^(-2) ), 1/d, nu)
    }
    
    # Compute quantiles located at (-Inf,0] and at (0,+Inf)
    F0 = pgat(0, mean = 0, sd = 1, nu = nu, d = d, xi = xi)
    n = length(p)
    result = rep(NA,n)
    indexLessThanF0 = which (p <= F0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanF0)
    if(sizeIndex == 0) {
      
        U = Up (p = p, nu = nu, d = d, xi = xi)
        result = ( U * nu / (1 - U) )^( 1/d ) * xi
        
    } else if (sizeIndex == n) {
      
         L = Lp (p = p, nu = nu, d = d, xi = xi)
        result = - ( nu/L - nu )^( 1/d ) * 1/xi
        
    } else if (TRUE) {
      
        L = Lp (p = p[indexLessThanF0], nu = nu, d = d, xi = xi)
        U = Up (p = p[-indexLessThanF0], nu = nu, d = d, xi = xi)
        result[indexLessThanF0] = - ( nu/L - nu )^( 1/d ) * 1/xi
        result[-indexLessThanF0] = ( U * nu / (1 - U) )^( 1/d ) * xi 
    }
    
    # Return Value:
    result * sd + mean
  }


#' @export
rgat <-  
  function(n, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)  
  {   
    
    randomUnif = runif(n = n, min = 0, max = 1)
    result = qgat(p = randomUnif, mean = mean, sd = sd, nu = nu, d = d, xi = xi)
    
    # Return Value:
    result
  }


gat.valid.pars = function(mean, sd, nu , d, xi){
  if(!(xi > 0) || !(d > 0) || !(nu > 0) || !(sd > 0) )
    return(FALSE)
  return(TRUE)
}


#' Estimate GAT parameters
#'
#' Functions to estimate all parameters of the GAT distribution from a vector
#' of iid observations. 
#'
#' It optimizes the log-likelihood based on dgat data
#' @name gat.fit
#' @rdname gat.fit
#' @param x Numeric vector of observations for estimating parameters.
#' @param start (optional) starting values of parameters as c(mean,sd,nu,d,xi) 
#' for the optimization algorithm 
#' @param lower.bound (optional) lower bounds for the optimization algorithm 
#' @param upper.bound (optional) lower bounds for the optimization algorithm 
#'
#' @return
#' An object with the optimization output, including estimated parameters,
#'   convergence status, objective value, and additional diagnostic information. The
#'   optimization is performed using \code{\link[Rsolnp]{solnp}}.
#' @author Thiago do Rego Sousa
#'
#' @examples
#' # simulate random values from GAT distribution
#' x = rgat(n = 1000, mean = 2, sd = 1, nu = 2, d = 1, xi = 3)
#' # estimate the parameters using the observations x
#' gat.fit(x)$pars
#'
#' @export
gat.fit <- function(x, start, lower.bound, upper.bound, control = NULL)  {   
  
    if(missing(start)){
      start = c(median(x),mad(x),1,1,1)
    }
    if(missing(lower.bound)){
      lower.bound = c(median(x) - 2*mad(x), 0.01, 0.01, 0.01, 0.01)
    }
    if(missing(upper.bound)){
      upper.bound = c(median(x) + 2*mad(x), 10, 10, 10, 10)
    }
    
    llh = function(pars){
      
      mean = pars[1]
      sd = pars[2]
      nu = pars[3]
      d  = pars[4]
      xi = pars[5]

      if(gat.valid.pars(mean,sd, nu , d, xi))   
        return(-sum(log(dgat(x = x, mean = mean, sd = sd, nu = nu, d = d, xi = xi)))) 
      else
        print('here')
        return(1e99)
    }
    
    fit <- solnp(pars = start, fun = llh, 
                  LB = lower.bound, UB = upper.bound, control = control)
    
    # Return Value:
    return(fit)
}

