\name{GEVStableGarch-package}


\alias{GEVStableGarch-package}
\alias{GEVStableGarch}


\docType{package}


\title{

ARMA-GARCH/APARCH modelling with GEV and stable distributions

}


\description{

This package is designed to perform maximum likelihood estimation of ARMA-GARCH/APARCH models 
with GEV, stable and generalized asymmetric t (GAt). Other common conditional distributions such as the normal, skew Student's t and GED are also allowed. 

}


\details{

    \tabular{ll}{
      Package: \tab GEVStableGarch \cr
      Type: \tab Package \cr
      Version: \tab 1.1 \cr
      Date: \tab 2015-07-19 \cr
      License: \tab GPL(>=2) \cr
      Depends: \tab R(>= 2.15.0), fGarch, fExtremes, stabledist, skewt, Rsolnp \cr
    }


GARCH models have proven to be highly effective for analyzing financial data over the past decades. In particular, the combination of ARMA-GARCH models with stable and generalized extreme distribution (GEV) distributions was successfully applied for forecasting volatility and for the measurement of Value at Risk (VAR).

Choosing the normal distribution as probability distribution for the innovations was a common choice in the beginning of the development of ARCH-type models. But recent research Nolan (1999), Mittnik et al. (2002), Mittnik and Paolella (2003), Curto et al. (2006), Frain (2009), Zhao et al. (2011) has shown that other distributions should be considered, specially because normal distribution can not account for fat tails and asymmetry found in real data.

This package contains functions for estimating and simulating combined ARMA-GARCH or ARMA-APARCH models with error distributions following GEV or stable densities. Our implementation also allows the estimation using maximum log-likelihood under different assumptions, Normal, Student's t, skew Student's t from Fernandez and Steel (1998), Generalized Asymmetric t (GAt) from Paolella (1997) and Generalized Error Distribution (GED).

The current version of package \pkg{GEVStableGarch} has a new algorithm that allows the user to enforce stationarity during estimation. Aditionally, it contains functions for selecting the best model according to a goodness-of-fit criteria.

}


\section{Time Series Simulation}{

  contains functions to simulate ARMA-GARCH/APARCH processes with conditional GEV or stable
  distributions. Note: These routines were adapted from /code{garchSpec} and /code{garchSim} functions available in \pkg{fGarch} package.
  
  Functions:
  
  \tabular{ll}{
    \code{gsSpec} \tab Specifies an univariate ARMA-GARCH/APARCH model, \cr 
    \code{gsSim} \tab Simulates a ARMA-GARCH/APARCH process. }
  
}


\section{Parameter Estimation}{

  contains functions to fit the parameters of ARMA-GARCH/APARCH
  time series processes.
  
  Functions:
  
  \tabular{ll}{
    \code{gsFit} \tab Fits the parameters of a GARCH process. This function also
    provides an algorithm to enforce stationarity during estimation, \cr 
    \code{gsSelect} \tab Selects the best model according to a goodness-of-fit
    criteria. } 
}


\section{Other Conditional Distribution Functions}{

  contains functions to compute density, distribution, quantiles and generate random values using important conditional distributions used in the garch literature.
  
  Functions:
  
  \tabular{ll}{
    \code{[dpqr]skstd} \tab Skew Student's t distribution function from Fernandez and Steel (1998), \cr 
    \code{[dpqr]gat} \tab Generalized Asymmetric t distribution (GAt) from Paolella (1997). Also known as t3-distribution, \cr 
    \code{gsMomentAparch} \tab Computes APARCH Moments for the following distributions: stable, GEV, GAt, normal, skew Student's t from Fernandez and Steel (1998) and GED.
  } 
  
}



\author{

Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano 
and Silvia Regina Costa Lopes.

Maintainer: Thiago do Rego Sousa <thiagoestatistico@gmail.com>

}


\keyword{ stable distribution }
\keyword{ GEV distribution }
\keyword{ GAt distribution }
\keyword{ t3-distribution }
\keyword{ ARMA-GARCH }
\keyword{ ARMA-APARCH }