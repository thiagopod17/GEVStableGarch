\name{GSgarch.Sim}
\alias{GSgarch.Sim}
\title{
Simulation of ARMA-GARCH/APARCH process with GEV and stable conditional distributions
}
\description{
This functions simulate time series following ARMA-GARCH/APARCH models with several conditional 
distributions, including GEV and stable distributions.
}
\usage{
GSgarch.Sim(N = 1000, mu = 0.1, a = c(0.5, 0.3), b = c(-0.4, 0.3, -0.1), 
omega = 0.05, alpha = c(0.1), gm = c(0), beta = c(0.1, 0.05, 0.03), 
delta = 2, skew = 0, shape = 3, cond.dist = "norm")
}
\arguments{
  \item{N}{
	The size of simulated time series
}
  \item{mu,a,b}{
	These are the parameters of the ARMA part of the model to be simulated.
}
\item{omega,alpha,beta}{
	GARCH parameters of the model to simulate. Both alpha and beta are given as vectors.
}
\item{gm,delta}{
	APARCH parameters of the model to simulate. Parameters gm and alpha must have the same size. 
}
  \item{skew}{
	The skew parameter of the conditional distribution. Used only when cond.dist equals \dQuote{stable} or \dQuote{sstd}.
}
  \item{shape}{
	Shape parameter of the conditional distribution. Used only when cond.dist equals to
	\dQuote{t-student}, \dQuote{skew t-student},\dQuote{GEV} or \dQuote{stable}.
}
  \item{cond.dist}{
	The conditional distribution of the model to be estimated. This parameter should be one of the following lowercase strings:
	\dQuote{norm} (Normal Distribution), \dQuote{t-student} (t-Student distribution), \dQuote{sstd} (skew t-Student Distribution), 
	\dQuote{GEV} (Generalized Extreme Value Distribution) or \dQuote{stable} (Stable Distribution). 
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
x <- GSgarch.Sim(N = 2500, mu = 0.1,a = c(0.2,0.3),b = c(0.2,0.5),
omega = 0.1, alpha = c(0.1,0.2),beta = c(0.1,0.1),gm=c(0.3,-0.3),
delta = 1,skew = 0.3,shape = 1.9, cond.dis = "stable")
}