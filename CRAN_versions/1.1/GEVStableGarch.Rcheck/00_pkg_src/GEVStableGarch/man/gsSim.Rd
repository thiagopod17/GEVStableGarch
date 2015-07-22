\name{gsSim}
\alias{gsSim}
\title{
Simulation of ARMA-GARCH/APARCH process
}
\description{
This functions simulate time series following ARMA-GARCH/APARCH models with several conditional distributions, including GEV and stable distributions.}
\usage{
gsSim(spec = garchSpec(), n = 100, n.start = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spec}{
A model specified with function gsSpec
}
  \item{n}{
The size of simulated time serie
}
  \item{n.start}{
%%     ~~Describe \code{n.start} here~~
}
}
\details{
The initial values of the time series
are fixed and the recursion formulas of the model are used to simulate the dynamics of the
process. We do not verify the stationarity conditions of the model because the simulation of
non-stationary process could also be of interest.
}
\value{
  The function returns an object containing the following items:
	\item{model }{A string describing the estimated model.}
	\item{cond.dist }{The conditional distribution used to fit the model.}
	\item{series}{An array of two columns. The first column is the simulated process \deqn{X_t}
	and the second one is the \deqn{\sigma_t} simulated process.}
}
\references{
Brockwell, P. J, e Davis, R. A. Introduction to Time Series and Forecasting.
Springer, New York, 1996.

Nolan, J. P. Numerical calculations of stable densities and distribution
functions. Communications in Statistics - Stochastic Models 13: 759-
774, 1997.
}
\author{
Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano 
and Silvia Regina Costa Lopes.
}
\examples{
# Simulation of a ARMA-APARCH process with stable conditional distribution
#x <- GSgarch.Sim(N = 2500, mu = 0.1,a = c(0.2,0.3),b = c(0.2,0.5),
#omega = 0.1, alpha = c(0.1,0.2),beta = c(0.1,0.1),gm=c(0.3,-0.3),
#delta = 1,skew = 0.3,shape = 1.9, cond.dis = "stable")
}