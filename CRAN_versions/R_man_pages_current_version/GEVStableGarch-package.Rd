\name{GEVStableGarch-package}
\alias{GEVStableGarch-package}
\alias{GEVStableGarch}
\docType{package}
\title{
ARMA-GARCH or ARMA-APARCH modelling with GEV and stable distributions
}
\description{
This package is design to perform maximum likelihood estimation of ARMA-GARCH/APARCH models 
with Generalized Extreme Distribution (GEV) or stable conditional distributions.
}
\details{
\tabular{ll}{
Package: \tab GEVStableGarch\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2015-07-19\cr
License: \tab GPL(>=2)\cr
Depends: \tab R(>= 2.15.0), fGarch, fExtremes, stabledist, skewt, Rsolnp\cr
}
This package contains functions for estimating and simulating combined ARMA-GARCH or ARMA-APARCH models with error distributions following GEV or stable densities. 
The current version has a new algorithm that allows the user to enforce stationarity during estimation. 
The package also allows the user to estimate
several models and choose the best one according to a goodness-of-fit Criteria.
}
\author{
Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano 
and Silvia Regina Costa Lopes

Maintainer: Thiago do Rego Sousa <thiagoestatistico@gmail.com>
}
\references{
Brockwell, P. J, e Davis, R. A. Introduction to Time Series and Forecasting.
Springer, New York, 1996.

Mittnik S. , Paolella M.S. , Rachev, S.T. Stationarity of stable power-
GARCH processes. Journal of Econometrics 106, 97-107. 2002.

Nolan, J. P. Numerical calculations of stable densities and distribution
functions. Communications in Statistics - Stochastic Models 13: 759-
774, 1997.

Nolan, J. P. Maximum likelihood estimation and diagnostics for stable distribution. 
In O. E.Barndorff-Nielsen, T. Mikosch, e S. I. Resnick (Eds.), 
Levy Processes: Theory and Applications: 379-400. Boston, Birkhauser, 2000.

Wurtz, D., Y. Chalabi, e L. Luksan. Parameter estimation of ARMA
models with GARCH/APARCH errors: An R and SPlus software implementation.
Journal of Statistical Software, 2006.

Zhao X. ,Scarrott C. J. , Oxley L. , Reale M. GARCH dependence in extreme
value models with Bayesian inference. Mathematics an Computers
in Simulation, Vol. 81, Issue 7, 1430-1440. 2011.
}
\keyword{ stable distribution }
\keyword{ GEV }
\keyword{ GAt }
\keyword{ t3-distribution }
\keyword{ ARMA-GARCH }
\keyword{ ARMA-APARCH }