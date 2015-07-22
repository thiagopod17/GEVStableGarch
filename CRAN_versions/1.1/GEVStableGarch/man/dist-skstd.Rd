\name{skstd}
\alias{skstd}

\alias{dskstd}
\alias{pskstd}
\alias{qskstd}
\alias{rskstd}


\title{
Skew Student's t Distribtuion from Fernandez and Steel (2002)
}


\description{

    Functions to compute density, distribution function, 
    quantile function and to generate random variates 
    for the Skew Student's t Distribtuion from Fernandez and Steel (2002) withouth reparemeterization. 
    
}
\usage{
dskstd(x, mean = 0, sd = 1, nu = 3, xi = 1, log = FALSE)
pskstd(q, mean = 0, sd = 1, nu = 3, xi = 1)
qskstd(p, mean = 0, sd = 1, nu = 3, xi = 1)
rskstd(n, mean = 0, sd = 1, nu = 3, xi = 1)
}

\arguments{
    \item{mean, sd, nu, xi}{
        location parameter \code{mean},
        scale parameter \code{sd},
        shape parameter \code{nu}.
        asymmetry parameter \code{xi}.
        }
        
    \item{n}{
        the number of observations.
        } 
    \item{p}{
        a numeric vector of probabilities.
        }
    \item{x, q}{
        a numeric vector of quantiles.
        }
    \item{log}{
        a logical; if TRUE, densities are given as log densities.
        }
}
\value{

    \code{d*} returns the density,
    \code{p*} returns the distribution function,
    \code{q*} returns the quantile function, and
    \code{r*} generates random deviates, \cr
    all values are numeric vectors.
 
}



\references{
%% ~put references to the literature/web site here ~
}
\author{
Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano 
and Silvia Regina Costa Lopes.
}

\examples{
%% asdf
}
\keyword{ distribution }