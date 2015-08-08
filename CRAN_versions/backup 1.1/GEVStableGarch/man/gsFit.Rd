\name{gsFit}
\alias{gsFit}


\title{

Estimation of ARMA-GARCH/APARCH models

}


\description{

This function uses Maximum Likelihood technique to estimate the parameters of
ARMA-GARCH or ARMA-APARCH model with several conditional distributions. The user can also enforce stationarity during estimation

}


\usage{

gsFit(formula = ~garch(1, 1), data, cond.dist = c("stableS0", "stableS1", 
"stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"), 
include.mean = TRUE, algorithm = c("sqp", "sqp.restriction", 
"nlminb", "nlminb+nm"), control = NULL, 
tolerance = NULL, title = NULL, description = NULL)

}


\arguments{

    \item{formula}{ 
        a formula object specifying mean and variance equation. The mean equation should be specified using the names: "arma"", "garch"" or "aparch". For example: ~arma(1,0)-garch(1,0) for AR(1)-ARCH(1), and ~arma(0,1)-garch(2,2) for MA(1)-GARCH(2,2). 
        }
        
   \item{data}{ 
        a numeric vector containing the data of
        the time series to be estimated.
        }

        \item{cond.dist}{ a character string naming conditional distribution of innovations. The package was created to accept the following distributions: \code{"stableS0"} (stable in S0-parameterization), \code{"stableS1"} (stable in S1-parameterization), \code{"stableS2"} (stable in S2-parameterization), \code{"gev"} and \code{"gat"}. Other common distributions are also possible such as the \code{"norm"}, \code{"std"}, \code{"sstd"} (skew Student's t defined in \pkg{fGarch} package), \code{"skstd"} (skew Student's t from Fernandez and Steel (1998)), \code{"ged"}.
        }
        
        \item{include.mean}{
          This is a boolean variable. It intercept is TRUE than we estimate the model with intercept.}
      
        \item{algorithm}{
      The algorithm to be used to search for the optimum value. The current version of the \pkg{GEVStableGarch} package implements
four different optimization procedures, namely the \code{"sqp"}, \code{"sqp.restriction"} (enforcing stationarity) , \code{"nlminb"} and \code{"nlminb+nm"}. See the details for more information about the implemented algorithms.
        
      }

        \item{control}{
        Control parameters used either in \dQuote{SQP} routine (see the RSolnp documentation) 
        or in \dQuote{nlminb} (see the nlminb R documentation).
        }

      \item{tolerance}{
        Tolerance for parameter estimation. When specified, the tolerance parameter should be passed as a list in the following format \code{list( TOLG = 1e-8, TOLSTABLE = 1e-2, TOLSTATIONARITY = 1e-3 )}, where \code{TOLG} is the tolerance used to set the boundary of the general model, \code{TOLSTABLE} is the boundary of parameters for stable distribution and \code{TOLSTATIONARITY} is the tolerance parameter used by the "sqp.restriction"" algorithm when searching for the stationary solution. 
        }

    \item{title}{ 
        a string with the title.
        }
        
    \item{description}{
        a string with a description.
        }
}


\details{

The starting values are key in getting any possible convergence value and they were chosen to reflect this. For example, the GEV shape starting value was chosen to be 0.01 since in many cases we can not obtain good results.
There is also the possibility to fit models with conditional
\dQuote{Normal}, \dQuote{t-student} and \dQuote{skew t-student} distributions but the main 
purpose of this routine is to deal with \dQuote{GEV} and \dQuote{stable}. 
The parameters will be interpreted according to the following 
equations (see Wurtz et al. ,2009)
  \deqn{				
				X_t =  \mu + \sum_{i=1}^m a_i X_{t-i} + \sum_{j=1}^n b_j \varepsilon_{t-j} + \varepsilon_t 
		}
	\deqn{	
				  \varepsilon_t = \sigma_t z_t,\;\;\;z_t  \stackrel{iid}{\sim} D(0,1)\;,\\
		}
	\deqn{	
				\sigma_t^\delta   = \omega + \sum_{i=1}^p \alpha_i(\varepsilon_{t-i}-\gamma_i|\varepsilon_{t-i}|)^\delta + \sum_{j=1}^q \beta_j \sigma_{t-j}^\delta		
		}
where \deqn{\mathcal{D}_{\vartheta}(0,1)} is the density of the innovations with zero location and unit scale and \deqn{\vartheta} are additional distributional parameters that describe the skew and the shape of the distribution.

ESTIMATION ALGORITHM
Most software packages implement the estimation of GARCH models without imposing sta- tionarity, but restricting the parameter set by appropriate bounds. This last approach was implemented in the GEVStableGarch package through the following algorithms: "sqp", "nlminb" and "nlminb + nm". The first two algorithms search for the optimum value by restricting the parameter set to appropriate lower and upper bounds. The last implements a two step optimization procedure, which consists in starting the search by using the con- strained routine nlminb and then performing another search using an unconstrained method (in our case the Nelder-Mead method implemented in the R base function optim). 


Most software packages implement the estimation of GARCH models without imposing stationarity, but restricting the parameter set by appropriate bounds. This last approach was implemented in the \pkg{GEVStableGarch} package through the following algorithms: \code{"sqp"}, \code{"nlminb"} and \code{"nlminb + nm"}. The first two algorithms search for the optimum value by restricting the parameter set to appropriate lower and upper bounds. The last implements a two step optimization procedure, which consists in starting the search by using the constrained routine \code{nlminb} and then performing another search using an unconstrained method (in our case the Nelder-Mead method implemented in the \code{R} base function \code{optim}). This approach was suggested by Wuertz et al. (2009) since in many cases it leads to an improved solution (in terms of the likelihood function of the data).

Finally, the \code{"sqp.restriction"} algorithm gives us the possibility to maximize the log-likelihood function in such a way that the estimated model is still stationary. What made this possible was the availability of the APARCH moments formulas throughout the literature for many conditional distributions and the derivation of the stable case given in this paper.

}


\value{
  
    \code{gsFit} \cr
    returns a S4 object of class \code{"GEVSTABLEGARCH"} with the following slots: 
  
    \item{@call}{
        the call of the \code{gsFit} function.
    }
    \item{@formula}{
        a list with two formula entries, one for the mean and the other 
        one for the variance equation.
    }
    \item{@method}{
        a string denoting the optimization method.
    }  
   \item{@data}{ 
        a numeric vector containing the data of the estimated time sereis.
        }
        
    \item{@fit}{ 
    a list with the results from the parameter estimation:
      \code{par} - the estimated parameters;
      \code{llh} - the estimated negative log-likelihood function;
      \code{hessian} - the hessian matrix returned by the optimization algorithm;
      \code{ics} - the value of the goodness-of-fit measures (AIC, AICc and BIC) (See Brockwell and Davis, 2002 for more details);
      \code{order} - a list with the ARMA and garch ORDERS;
      \code{cond.dis} - the conditional distribution;
      \code{se.coef} - standard errors of the estimated parameters;
      \code{tValue} - tValue of the estimated parameters;
      \code{matcoef} - an organized matrix with the estimated parameters.
    }
    
    \item{@residuals}{ 
    a numeric vector with the residual values (\eqn{\varepsilon_t}).
    }

    \item{@h.t}{
        a numeric vector with the residual values
        (\eqn{h_t = \sigma_t^\delta}). 
    }
    
    \item{@sigma.t}{
        a numeric vector with the conditional standard deviation \eqn{\sigma_t}.
    }
    
    \item{@title}{
        a string with the title.
    }
    
    \item{@description}{
        a string with a description.
    } 
  
    The entries of the @fit slot show the results from the
    optimization.

}


\section{Enhances}{

The estimation of ARMA-GARCH/APARCH models with conditional stable distribution is mainly dependent
on the time taken during the calculation of density points. Although the package can be used with the R \pkg{stabledist} package, the estimation of such models is only feasible if we use the fast implementation of stable densities using the R \pkg{stable} package designed by Nolan. The \pkg{stable} package implements a faster computation of stable 
densities that are accurately enough to perform numerical optimization.
Package \pkg{stable} is available at \url{http://www.robustanalysis.com}.

}


\references{

    Brockwell, P.J., Davis, R.A. (1996).
    \emph{Introduction to Time Series and Forecasting}.
    Springer, New York.

    Mittnik, S., Paolella, M.S., Rachev, S.T. (2002).
    \emph{Stationarity of stable power-GARCH processes}.
    Journal of Econometrics 106, 97--107.

    Nolan, J.P. (1997).
    \emph{Numerical calculations of stable densities and distribution functions}.
    Communications in Statistics - Stochastic Models 13, 759--774.

    Nolan, J.P. (1997). 
    \emph{Maximum likelihood estimation and diagnostics for stable distribution}. 
    In O. E.Barndorff-Nielsen, T. Mikosch, e S. I. Resnick (Eds.),
    Levy Processes: Theory andApplications,
    Boston, Birkhauser, 379--400.

    Wuertz, D., Chalabi, Y., with contribution from Miklovic, M., Boudt, C., Chausse, P., and others (2013).
    \emph{fGarch: Rmetrics - Autoregressive Conditional Heteroskedastic Modelling, R package version 3010.82},   
    \url{http://CRAN.R-project.org/package=fGarch}.

    Wuertz, D., Chalabi, Y., Luksan, L. (2009).
    \emph{Parameter Estimation of ARMA Models with GARCH/ APARCH Errors: An R and SPlus SoftwareImplementation}.   
    Journal of Statistical Software, forthcoming,
    \href{http://www-stat.wharton.upenn.edu/~steele/Courses/956/RResources/GarchAndR/WurtzEtAlGarch.pdf}{http://www-stat.wharton.upenn.edu/~steele/...WurtzEtAlGarch.pdf}.

    Zhao, X., Scarrott, C.J., Oxley, L., Reale, M. (2011).
    \emph{GARCH dependence in extreme value models with Bayesian inference}.
    Mathematics an Computers in Simulation, 81, Issue 7, 1430--1440.

    Yinyu, Ye. (1987).
    \emph{Interior Algorithms for Linear, Quadratic, and Linearly Constrained Non-Linear Programming}.
    Ph.D. Thesis, Department of EES,
    Stanford University.

}

\author{

Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano
and Silvia Regina Costa Lopes

}


\seealso{

\code{\link{gsSelect}}

}

\examples{
# This examples uses the dataset of the package fGarch to estimate
# an ARMA(1,1)-GARCH(1,1) with GEV conditional distribution.
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]
gev.model = gsFit(data = x , formula = ~garch(1,1), cond.dist = "norm")
}



