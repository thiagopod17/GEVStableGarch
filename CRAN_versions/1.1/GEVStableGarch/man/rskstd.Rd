\name{rskstd}
\alias{rskstd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rskstd(n, mean = 0, sd = 1, nu = 3, xi = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (n, mean = 0, sd = 1, nu = 3, xi = 1) 
{
    randomUnif = runif(n = n, min = 0, max = 1)
    result = qskstd(p = randomUnif, mean = mean, sd = sd, nu = nu, 
        xi = xi)
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
