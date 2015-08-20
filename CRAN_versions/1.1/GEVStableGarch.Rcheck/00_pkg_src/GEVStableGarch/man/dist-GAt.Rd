\name{gat}
\alias{gat}
\alias{dgat}
\alias{pgat}
\alias{qgat}
\alias{rgat}


\title{

Generalized Asymmetric t Distribution

}


\description{

    Functions to compute density, distribution function, 
    quantile function and to generate random variates 
    for the generalized asymmetric t distribution (GAt). Notice that
    this is the same t3-distribution mentioned in the literature and defined 
    by Paolella (1997) and Mittnik and Paolella (2000). 
    The GAt distribution includes the Student's t, Laplace, 
    Cauchy and the normal distribution when the shape parameter (\eqn{\nu \rightarrow \infty}) (see Mittnik and Paolella (2000)).
    
}


\usage{

dgat(x, mean = 0, sd = 1, nu = 2, d = 3, xi = 1, log = FALSE)
pgat(q, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
qgat(p, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
rgat(n, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)

}


%- maybe also 'usage' for other objects documented here.
\arguments{

    \item{mean, sd, d, nu, xi}{
        location parameter \code{mean},
        scale parameter \code{sd},
        shape 1 parameter \code{nu},
        shape 2 parameter \code{d},
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

    Mittnik, S., Paolella, M.S. (2000).
   \emph{Prediction of Financial Downside-Risk with Heavy-Tailed Conditional Distributions}
    Available at SSRN 391261.

    Paolella, M. (1997).
    \emph{Tail Estimation and Conditional Modeling of Heteroskedstic Time-Series}.
    Ph.D Thesis, Institute of Statistics and Econometrics,
    Christian Albrechts University of Kiel.

}


\author{

Thiago do Rego Sousa

}


\examples{

# Simulate Random Values and compare with
# the empirical density and probability functions
# Note: This example was addapted from "sstd {fGarch} R Documentation"

# Configure plot and generate random values
par(mfrow = c(2, 2))
set.seed(1000)
r = rgat(n = 1000)
plot(r, type = "l", main = "GAt Random Values", col = "steelblue")

# Plot empirical density and compare with true density:
hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue")
box()
x = seq(min(r), max(r), length = 201)
lines(x, dgat(x), lwd = 2)

# Plot density function and compare with true df:
plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue",
     ylab = "Probability")
lines(x, pgat(x), lwd = 2)

# Compute quantiles:
# Here we compute the quantiles corresponding to the probability points from 
# -10 to 10 and expect to obtain the same input sequence
round(qgat(pgat(q = seq(-10, 10, by = 0.5))), digits = 6)

}


\keyword{ distribution }
