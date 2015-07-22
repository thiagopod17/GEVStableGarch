\name{gsFit}
\alias{gsFit}

\title{
Estimation of ARMA-GARCH/APARCH models
}
\description{
This function uses Maximum Likelihood technique to estimate the parameters of
ARMA-GARCH or ARMA-APARCH model with several conditional distributions. The user can also enforce stationarity during estimation}
\usage{
gsFit(formula = ~garch(1, 1), data, cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "GAt", "norm", "std", "sstd", "skstd", "ged"), include.mean = TRUE, algorithm = c("sqp", "sqp.restriction", "nlminb", "nlminb+nm"), control = NULL, tolerance = NULL, title = NULL, description = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{cond.dist}{
The conditional distribution of the model to be estimated. This parameter should be one of the following strings:
  \dQuote{norm} (Normal Distribution), \dQuote{t-student} (t-Student distribution), \dQuote{sstd} (skew t-Student Distribution), 
	\dQuote{GEV} (Generalized Extreme Value Distribution) or \dQuote{stable} (Stable Distribution)
}
  \item{include.mean}{
This is a boolean variable. It intercept is TRUE than we estimate the model with intercept, otherwise we will not use the intercept variable.}
  \item{algorithm}{
The algorithm to be used to optimize the likelihood of the process. This parameter should be one of the 
  following strings: \dQuote{nlminb} ( nlminb R internal Routine for optimization) or \dQuote{SQP} 
	( Sequential quadratic programming routine implemented in package RSolnp).
}
  \item{control}{
Control parameters used either in \dQuote{SQP} routine (see the RSolnp documentation) 
or in \dQuote{nlminb} (see the nlminb R documentation).
}
  \item{tolerance}{
Tolerance for parameter estimation.
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{description}{
%%     ~~Describe \code{description} here~~
}
}
\details{
There is also the possibility to fit models with conditional
\dQuote{Normal}, \dQuote{t-student} and \dQuote{skew t-student} distributions but the main 
purpose of this routine is to deal with \dQuote{GEV} and \dQuote{stable}. 
The parameters will be interpreted according to the following 
equations (see Wurtz et al. ,2006)
  \deqn{				
				X_t =  \mu + \sum_{i=1}^m a_i X_{t-i} + \sum_{j=1}^n b_j \varepsilon_{t-j} + \varepsilon_t 
		}
	\deqn{	
				  \varepsilon_t = \sigma_t z_t,\;\;\;z_t  \stackrel{iid}{\sim} D(0,1)\;,\\
		}
	\deqn{	
				\sigma_t^\delta   = \omega + \sum_{i=1}^p \alpha_i(\varepsilon_{t-i}-\gamma_i|\varepsilon_{t-i}|)^\delta + \sum_{j=1}^q \beta_j \sigma_{t-j}^\delta		
		}
where \deqn{X_t} is the ARMA process, \deqn{\sigma_t} is the APARCH process and \deqn{D(0,1)} 
is the standard density distribution to be used in the model.
}
\value{
	The function returns an object containing the following items:
	\item{order }{The estimated parameter set for the chosen model.}
	\item{hessian }{The estimated Hessian matrix.}
	\item{model }{A string describing the estimated model.}
	\item{cond.dist }{The conditional distribution used to fit the model.}
	\item{data }{The time series data used to fit the model.}     
	\item{llh }{The negative log likelihood of the estimated model}
	\item{par }{a vector containing the estimated parameters}     
    \item{hessian }{The estimated Hessian matrix.}  
	\item{ARMA.res }{ARMA residuals }  
	\item{GARCH.sig }{GARCH or APARCH residuals.} 
	\item{aic }{The reported AIC (Akaike Information Criterion) for the current model.}      
	\item{aicc }{The reported corrected AIC (AICc) for the current model.}      
	\item{bic }{The reported BIC (Bayesian Information Criterion. }      
	\item{se.coef }{Standard error for the estimated parameters.}  
	\item{tval}{Calculated t-value for the estimated parameters.}     
	\item{matcoef}{An organized matrix to present the objects.}  
}
\section{Enhances}{
The estimation of ARMA-GARCH/APARCH models with conditional stable distribution is mainly dependent
on the time taken during the calculation of density points. Our routine uses the standard R stable density
implemented in package \pkg{stabledist}. There is also another numerical routine to calculate stable 
densities which is part of an R package called \pkg{stable}. This package implements a faster computation of stable 
densities that are accurately enough to perform numerical optimization. 
When \pkg{GEVStableGarch} loads it verifies if package \pkg{stable} is available. 
If it is found the function GSGarch.Fit will perform optimization using \pkg{stable}'s routine. 
Otherwise it will use stable density from package \pkg{stabledist}.
Package \pkg{stable} is available at \samp{http://www.robustanalysis.com}.
}
\references{
Mittnik S. , Paolella M.S. , Rachev, S.T. Stationarity of stable power-
GARCH processes. Journal of Econometrics 106, 97-107. 2002.

Wurtz, D., Y. Chalabi, e L. Luksan. Parameter estimation of ARMA
models with GARCH/APARCH errors: An R and SPlus software implementation.
Journal of Statistical Software, 2006.

Zhao X. ,Scarrott C. J. , Oxley L. , Reale M. GARCH dependence in extreme
value models with Bayesian inference. Mathematics an Computers
in Simulation, Vol. 81, Issue 7, 1430-1440. 2011.

Yinyu Ye. Interior Algorithms for Linear, Quadratic, and
Linearly Constrained Non-Linear Programming. Ph.D. Thesis, Department
of EES, Stanford University, 1987.
}
\author{
Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano 
and Silvia Regina Costa Lopes
}
\examples{
# This examples uses the dataset of the package fGarch to estimate
# an ARMA(1,1)-GARCH(1,1) with GEV conditional distribution.
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]
gev.model = gsFit(data = x , formula = ~garch(1,1), cond.dist = "norm")
}
