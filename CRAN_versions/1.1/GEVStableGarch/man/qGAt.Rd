\name{qGAt}
\alias{qGAt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qGAt(p, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{mean}{
%%     ~~Describe \code{mean} here~~
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{xi}{
%%     ~~Describe \code{xi} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (p, mean = 0, sd = 1, nu = 2, d = 3, xi = 1) 
{
    Lp <- function(p = p, nu = nu, d = d, xi = xi) {
        qbeta((1 + xi^2) * p, nu, 1/d)
    }
    Up <- function(p = p, nu = nu, d = d, xi = xi) {
        qbeta((p - 1/(1 + xi^2)) * (1 + xi^(-2)), 1/d, nu)
    }
    F0 = pGAt(0, mean = 0, sd = 1, nu = nu, d = d, xi = xi)
    n = length(p)
    result = rep(NA, n)
    indexLessThanF0 = which(p <= F0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanF0)
    if (sizeIndex == 0) {
        U = Up(p = p, nu = nu, d = d, xi = xi)
        result = (U * nu/(1 - U))^(1/d) * xi
    }
    else if (sizeIndex == n) {
        L = Lp(p = p, nu = nu, d = d, xi = xi)
        result = -(nu/L - nu)^(1/d) * 1/xi
    }
    else if (TRUE) {
        L = Lp(p = p[indexLessThanF0], nu = nu, d = d, xi = xi)
        U = Up(p = p[-indexLessThanF0], nu = nu, d = d, xi = xi)
        result[indexLessThanF0] = -(nu/L - nu)^(1/d) * 1/xi
        result[-indexLessThanF0] = (U * nu/(1 - U))^(1/d) * xi
    }
    result * sd + mean
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
