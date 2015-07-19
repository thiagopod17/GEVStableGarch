\name{pGAt}
\alias{pGAt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pGAt(x, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, mean = 0, sd = 1, nu = 2, d = 3, xi = 1) 
{
    if (length(mean) == 5) {
        xi = mean[5]
        d = mean[4]
        nu = mean[3]
        sd = mean[2]
        mean = mean[1]
    }
    if (sd <= 0 || nu <= 0 || xi <= 0 || d <= 0) 
        stop("Failed to verify condition:\n           sd <= 0 || nu <= 0 || xi <= 0 || d <= 0")
    L <- function(z, nu = nu, d = d, xi = xi) {
        nu/(nu + (-z * xi)^d)
    }
    U <- function(z, nu = nu, d = d, xi = xi) {
        pw = (z/xi)^d
        return(replace(pw/(nu + pw), which(pw == Inf, arr.ind = TRUE), 
            1))
    }
    z = (x - mean)/sd
    n = length(z)
    arg = z
    indexLessThanZero = which(z <= 0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanZero)
    if (sizeIndex == 0) {
        arg = 1/(1 + xi^2) + 1/(1 + xi^(-2)) * pbeta(U(z = arg, 
            nu = nu, d = d, xi = xi), 1/d, nu)
    }
    else if (sizeIndex == n) {
        arg = 1/(1 + xi^2) * pbeta(L(z = arg, nu = nu, d = d, 
            xi = xi), nu, 1/d)
    }
    else if (TRUE) {
        arg[indexLessThanZero] = 1/(1 + xi^2) * pbeta(L(z = arg[indexLessThanZero], 
            nu = nu, d = d, xi = xi), nu, 1/d)
        arg[-indexLessThanZero] = 1/(1 + xi^2) + 1/(1 + xi^(-2)) * 
            pbeta(U(z = arg[-indexLessThanZero], nu = nu, d = d, 
                xi = xi), 1/d, nu)
    }
    arg
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
