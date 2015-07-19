\name{GSgarch.FitAIC}
\alias{GSgarch.FitAIC}
\title{
ARMA-GARCH/APARCH model with minimum AIC
}
\description{
This function estimates ARMA-GARCH/APARCH models with varying order and return the one with minimum AIC.
}
\usage{
GSgarch.FitAIC(data, mMAX = 1, nMAX = 1, pMAX = 1, qMAX = 1, 
cond.dist = "norm", algorithm = "sqp", APARCH = FALSE, 
intercept = TRUE, control = NULL)
}
\arguments{
   \item{data}{
	Data with the time series to be estimated. This object can not contain NULL elements.
}
  \item{mMAX}{
	Maximum m order of the ARMA(m,n) equation to be searched.
}
  \item{nMAX}{
	Maximum n  of the ARMA(m,n) equation order to be searched.
}
  \item{pMAX}{
	Maximum p order of the GARCH(p,q) equation to be searched.
}
  \item{qMAX}{
	Maximum q order of the GARCH(p,q) equation to be searched.
}
  \item{cond.dist}{
	The conditional distribution of the model to be estimated. This parameter should be one of the following strings:
	\dQuote{norm} (Normal Distribution), \dQuote{t-student} (t-Student distribution), \dQuote{sstd} (skew t-Student Distribution), 
	\dQuote{GEV} (Generalized Extreme Value Distribution) or \dQuote{stable} (Stable Distribution). 
}
  \item{algorithm}{
	The algorithm to be used to optimize the likelihood of the process. This parameter should be one of the 
	following strings: \dQuote{nlminb} ( nlminb R internal Routine for optimization) or \dQuote{SQP} 
	( Sequential quadratic programming routine implemented in package RSolnp).
}
  \item{APARCH}{
	You can set this boolean value as FALSE if we want to estimate a pure GARCH model or as TRUE if
	we want an APARCH model. 
}
  \item{intercept}{
	This is a boolean variable. It intercept is TRUE than we estimate the model with intercept, otherwise we
	will not use the intercept variable.
}
  \item{control}{
	Control parameters used either in \dQuote{SQP} routine (see the RSolnp documentation) 
	or in \dQuote{nlminb} (see the nlminb R documentation).
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
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]
GSgarch.FitAIC(data = x,1,0,1,0,cond.dist = "gev")
}