\name{skstd}
\alias{skstd}
\alias{dskstd}
\alias{pskstd}
\alias{qskstd}
\alias{rskstd}


\title{
Skew Student's t Distribtuion from Fernandez and Steel (1997)
}


\description{

    Functions to compute density, distribution function, 
    quantile function and to generate random values 
    for the Skew Student's t distribtuion from Fernandez and Steel (1997). Notice that this function is different from the skew Student's t \code{sstd} from package \pkg{fGarch}. Although the two distributions use the same approach from Fernandez and Steel (1997), the second one was reparameterized in such a way that when the parameters mean = 0 and sd = 1 the distribution will have a zero mean and unit variance, no matter the value of the skew parameter.
    The distributions \code{skstd} and \code{sstd} are the same when the asymmetry parameter \code{xi} equals to 1, in which case they correspond to the well known Student's t distribution \code{std}.
    
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
        shape parameter \code{nu},
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
    \code{r*} generates random deviates. \cr
    All values are numeric vectors.
 
}


\references{

    Fernandez, C., Steel, M.F. (1998).
    \emph{On Bayesian Modeling of Fat Tails and Skewness}.
    Journal of the American Statistical Association, Taylor & Francis Group, 93(441), 359--371.

}

\author{

Thiago do Rego Sousa.

}

\examples{

# Simulate Random Values and compare with
# the empirical density and probability functions
# Note: This example was addapted from "sstd {fGarch} R Documentation"

# Configure plot and generate random values
par(mfrow = c(2, 2))
set.seed(1000)
r = rskstd(n = 1000)
plot(r, type = "l", main = "Skew Student's t Random Values", col = "steelblue")

# Plot empirical density and compare with true density:
hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue")
box()
x = seq(min(r), max(r), length = 201)
lines(x, dskstd(x), lwd = 2)

# Plot density function and compare with true df:
plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue",
     ylab = "Probability")
lines(x, pskstd(x), lwd = 2)

# Compute quantiles:
# Here we compute the quantiles corresponding to the probability points from 
# -10 to 10 and expect to obtain the same input sequence
round(qskstd(pskstd(q = seq(-10, 10, by = 0.5))), digits = 6)

}
\keyword{ distribution }