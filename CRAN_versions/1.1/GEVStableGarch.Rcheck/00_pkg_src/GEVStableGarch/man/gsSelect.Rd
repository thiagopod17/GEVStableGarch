\name{gsSelect}
\alias{gsSelect}

\title{
Selects the best model according to goodness-of-fit
}
\description{
This function estimates ARMA-GARCH/APARCH models with varying order and return the one with minimum AIC, BIC or AICc
}
\usage{
gsSelect(data, order.max = c(1, 1, 1, 1), selection.criteria = c("AIC", "AICc", "BIC"), is.aparch = FALSE, cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "GAt", "norm", "std", "sstd", "skstd", "ged"), include.mean = TRUE, algorithm = c("sqp", "sqp.restriction", "nlminb", "nlminb+nm"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Data with the time series to be estimated. This object can not contain NULL elements.
}
  \item{order.max}{
Maximum order of models to search
}
  \item{selection.criteria}{
The goodness-of-fit measure to use be minimized.
}
  \item{is.aparch}{
Boolean variable indicating whether to search for ARMA-GARCH or ARMA-APARCH models
}
  \item{cond.dist}{
  The conditional distribution of the model to be estimated. This parameter should be one of the following strings:
	\dQuote{norm} (Normal Distribution), \dQuote{t-student} (t-Student distribution), \dQuote{sstd} (skew t-Student Distribution), 
	\dQuote{GEV} (Generalized Extreme Value Distribution) or \dQuote{stable} (Stable Distribution). 
}
  \item{include.mean}{
This is a boolean variable. It intercept is TRUE than we estimate the model with intercept, otherwise we
  will not use the intercept variable.
}
  \item{algorithm}{
  The algorithm to be used to optimize the likelihood of the process. This parameter should be one of the 
	following strings: \dQuote{nlminb} ( nlminb R internal Routine for optimization) or \dQuote{SQP} 
	( Sequential quadratic programming routine implemented in package RSolnp).
}
  \item{\dots}{
Additional parameters to be passed to function GSfit during estimation.
}
}
\value{
  The model with minimum AIC found within the specified range in parameters mMAX, nMAX, pMAX, qMAX. 
}
\references{
Brockwell, P. J, e Davis, R. A. Introduction to Time Series and Forecasting.
Springer, New York, 1996.
}
\author{
Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano 
and Silvia Regina Costa Lopes
}
\examples{
# AIC fit using models from ARMA(0,0)-GARCH(1,0) to ARMA(1,1)-GARCH(1,1)
# with GEV conditional distribution
#library(fGarch)
#data(dem2gbp)
#x = dem2gbp[, 1]
# GSgarch.FitAIC(data = x,1,0,1,0,cond.dist = "gev")
}