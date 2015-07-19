\name{garch.stationarity}
\alias{garch.stationarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
garch.stationarity(parm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
%%     ~~Describe \code{parm} here~~
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
function (parm) 
{
    omega <- parm[1 + m + n + 1]
    alpha <- parm[(2 + m + n + 1):(3 + m + n + p - 1)]
    gm <- parm[(2 + m + n + p + 1):(3 + m + n + p + p - 1)]
    beta <- parm[(2 + m + n + 2 * p + 1):(3 + m + n + 2 * p + 
        q - 1)]
    delta <- parm[2 + m + n + 2 * p + q + 1]
    skew <- parm[3 + m + n + 2 * p + q + 1]
    shape <- parm[(4 + m + n + 2 * p + q + 1):(4 + m + n + 2 * 
        p + q + lengthShape)]
    model = list(omega, alpha, beta, delta, skew, shape)
    .stationarityAparch(model = model, formula = formula, cond.dist = cond.dist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
