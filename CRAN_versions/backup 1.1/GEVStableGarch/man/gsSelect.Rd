\name{gsSelect}
\alias{gsSelect}


\title{

Selects the best model according to goodness-of-fit criteria

}


\description{

This function estimates ARMA-GARCH/APARCH models with varying order and returns the smallest goodness-of-fit criteria

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
      Maximum order of models to search. It must by a vector of the type
      \code{c(mMax, nMax, pMax, qMax)} with each entry representing the maximum order of the model that will be fitted when searching for the best model.
        }
        
      \item{selection.criteria}{
      The goodness-of-fit measure to use be minimized. Three different criterias are allowed: AIC, AICc or BIC. For more information see Brockwell and Davis (1996)
        }
        
      \item{is.aparch}{
      Boolean variable indicating whether to search for ARMA-GARCH or ARMA-APARCH models
        }
        
\item{cond.dist}{ a character string naming conditional distribution of innovations. The package was created to accept the following distributions: \code{"stableS0"}, \code{"stableS1"}, \code{"stableS2"}, \code{"gev"} and \code{"GAt"}. Other common distributions are also possible such as the \code{"norm"}, \code{"std"}, \code{"sstd"}, \code{"skstd"}, \code{"ged"}.
        }
        
      \item{include.mean}{
    This is a boolean variable. It intercept is TRUE than we estimate the model with intercept.
        }
        
        \item{algorithm}{
      The algorithm to be used to search for the optimum value. The current version of the \pkg{GEVStableGarch} package implements
four different optimization procedures, namely the \code{"sqp"}, \code{"sqp.restriction"} (enforcing stationarity) , \code{"nlminb"} and \code{"nlminb+nm"}. See the details for more information about the implemented algorithms.
        
      }
        
      \item{\dots}{
      Additional parameters to be passed to function \code{gsFit} during estimation.
        }
}


\value{
    Returns a S4 object of class "GEVSTABLEGARCH" with the best model. See the Value section of \code{gsFit} function for details. 
}


\references{

    Brockwell, P.J., Davis, R.A. (1996).
    \emph{Introduction to Time Series and Forecasting}.
    Springer, New York.

}


\author{

Thiago do Rego Sousa.

}


\examples{

# Best ARMA-GARCH model within the range ARMA(0,0)-GARCH(1,0) to ARMA(0,0)-GARCH(1,1)
# using the Corrected Akaike Information Criteria (AICc)
library(fGarch)
data(dem2gbp)
x = dem2gbp[,1]

model = gsSelect (data = x, order.max = c(0,0,1,1), is.aparch = FALSE, 
          algorithm = "sqp", cond.dist = "gev", selection.criteria = "AIC")

}