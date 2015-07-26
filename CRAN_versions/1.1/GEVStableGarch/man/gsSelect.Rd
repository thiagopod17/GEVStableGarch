\name{gsSelect}
\alias{gsSelect}


\title{

Selects the best model according to goodness-of-fit

}


\description{

This function estimates ARMA-GARCH/APARCH models with varying order and return the one with minimum AIC, BIC or AICc

}


\usage{
gsSelect(data, order.max = c(1, 1, 1, 1), 
selection.criteria = c("AIC", "AICc", "BIC"), is.aparch = FALSE, 
cond.dist = c("stableS0", "stableS1", "stableS2", "gev", 
"GAt", "norm", "std", "sstd", "skstd", "ged"), 
include.mean = TRUE, 
algorithm = c("sqp", "sqp.restriction", "nlminb", "nlminb+nm"), ...)
}


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
        
\item{cond.dist}{ a character string naming conditional distribution of innovations. The package was created to accept the following distributions: \code{"stableS0"}, \code{"stableS1"}, \code{"stableS2"}, \code{"gev"} and \code{"GAt"}. Other common distributions are also possible such as the \code{"norm"}, \code{"std"}, \code{"sstd"}, \code{"skstd"}, \code{"ged"}, but they are not the main feature of this work.
        }
        
      \item{include.mean}{
    This is a boolean variable. It intercept is TRUE than we estimate the model with intercept, otherwise we will not use the intercept variable.
        }
        
      \item{algorithm}{
      The algorithm to be used to optimize the likelihood of the process. This parameter should be one of the following strings: \dQuote{nlminb} ( nlminb R internal Routine for optimization) or \dQuote{SQP} 
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
library(fGarch)
data(dem2gbp)
x = dem2gbp[,1]

model = gsSelect (data = x, order.max = c(0,0,1,1), is.aparch = FALSE, 
          algorithm = "sqp", cond.dist = "norm", selection.criteria = "BIC")

}