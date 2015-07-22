\name{GAt}
\alias{GAt}
\alias{dGAt}
\alias{pGAt}
\alias{qGAt}
\alias{rGAt}

\title{
Generalized Asymmetric t Distribution
}
\description{

    Functions to compute density, distribution function, 
    quantile function and to generate random variates 
    for the generalized asymmetric t distribution. 
    
}


\usage{
dGAt(x, mean = 0, sd = 1, d = 2, xi = 1, log = FALSE)
pGAt(q, mean = 0, sd = 1, d = 2, xi = 1)
qGAt(p, mean = 0, sd = 1, d = 2, xi = 1)
rGAt(n, mean = 0, sd = 1, d = 2, xi = 1)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{mean, sd, d, xi}{
        location parameter \code{mean},
        scale parameter \code{sd},
        shape 1 parameter \code{d}.
        shape 2 parameter \code{xi}.
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
