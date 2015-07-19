\name{dGAt}
\alias{dGAt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dGAt(x, mean = 0, sd = 1, nu = 3, d = 2, xi = 1, log = FALSE)
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
  \item{log}{
%%     ~~Describe \code{log} here~~
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
function (x, mean = 0, sd = 1, nu = 3, d = 2, xi = 1, log = FALSE) 
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
    z = (x - mean)/sd
    n = length(z)
    arg = z
    indexLessThanZero = which(z < 0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanZero)
    if (sizeIndex == 0) {
        arg = arg/xi
    }
    else if (sizeIndex == n) {
        arg = -arg * xi
    }
    else if (TRUE) {
        arg[indexLessThanZero] = -arg[indexLessThanZero] * xi
        arg[-indexLessThanZero] = arg[-indexLessThanZero]/xi
    }
    k = ((xi + 1/xi) * 1/d * nu^(1/d) * beta(1/d, nu))^(-1)
    result = (k * (1 + (arg^d)/nu)^(-nu - 1/d))/sd
    if (log) 
        result = log(result)
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
