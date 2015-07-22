\name{gsMomentAparch}
\alias{gsMomentAparch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsMomentAparch(cond.dist = c("stableS1", "gev", "GAt", "norm", "std", "sstd", "skstd", "ged"), shape = 1.5, skew = 0, delta = 1, gm = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cond.dist}{
%%     ~~Describe \code{cond.dist} here~~
}
  \item{shape}{
%%     ~~Describe \code{shape} here~~
}
  \item{skew}{
%%     ~~Describe \code{skew} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{gm}{
%%     ~~Describe \code{gm} here~~
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
function (cond.dist = c("stableS1", "gev", "GAt", "norm", "std", 
    "sstd", "skstd", "ged"), shape = 1.5, skew = 0, delta = 1, 
    gm = 0) 
{
    cond.dist = match.arg(cond.dist)
    if (cond.dist == "stableS1") 
        kappa = .stableS1MomentAparch(shape = shape, skew = skew, 
            delta = delta, gm = gm)
    if (cond.dist == "gev") 
        kappa = .gevMomentAparch(shape = shape, delta = delta, 
            gm = gm)
    if (cond.dist == "GAt") 
        kappa = .GAtMomentAparch(shape = shape, delta = delta, 
            skew = skew, gm = gm)
    if (cond.dist == "norm") 
        kappa = .normMomentAparch(delta = delta, gm = gm)
    if (cond.dist == "std") 
        kappa = .stdMomentAparch(shape = shape, delta = delta, 
            gm = gm)
    if (cond.dist == "sstd") 
        kappa = .sstdMomentAparch(shape = shape, skew = skew, 
            delta = delta, gm = gm)
    if (cond.dist == "skstd") 
        kappa = .skstdMomentAparch(shape = shape, skew = skew, 
            delta = delta, gm = gm)
    if (cond.dist == "ged") 
        kappa = .gedMomentAparch(shape = shape, delta = delta, 
            gm = gm)
    kappa
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
