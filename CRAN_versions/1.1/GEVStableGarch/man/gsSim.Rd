\name{gsSim}
\alias{gsSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsSim(spec = garchSpec(), n = 100, n.start = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spec}{
%%     ~~Describe \code{spec} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{n.start}{
%%     ~~Describe \code{n.start} here~~
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
function (spec = garchSpec(), n = 100, n.start = 100) 
{
    if (n < 2) 
        stop("The parameter 'n' must be > 2")
    stopifnot(class(spec) == "fGEVSTABLEGARCHSPEC")
    model = spec@model
    if (spec@rseed != 0) 
        set.seed(spec@rseed)
    n = n + n.start
    if (spec@distribution == "stableS0") 
        z = stabledist::rstable(n = n, alpha = model$shape, beta = model$skew, 
            pm = 0)
    if (spec@distribution == "stableS1") 
        z = stabledist::rstable(n = n, alpha = model$shape, beta = model$skew, 
            pm = 1)
    if (spec@distribution == "stableS2") 
        z = stabledist::rstable(n = n, alpha = model$shape, beta = model$skew, 
            pm = 2)
    if (spec@distribution == "gev") 
        z = rgev(n, xi = model$shape)
    if (spec@distribution == "GAt") 
        z = rGAt(n, nu = model$shape[1], d = model$shape[2], 
            xi = model$skew)
    if (spec@distribution == "norm") 
        z = rnorm(n)
    if (spec@distribution == "std") 
        z = rstd(n, nu = model$shape)
    if (spec@distribution == "sstd") 
        z = rsstd(n, nu = model$shape, xi = model$skew)
    if (spec@distribution == "skstd") 
        z = rskstd(n, nu = model$shape, xi = model$skew)
    if (spec@distribution == "ged") 
        z = rged(n, nu = model$shape)
    delta = model$delta
    z = c(rev(spec@presample[, 1]), z)
    h = c(rev(spec@presample[, 2]), rep(NA, times = n))
    y = c(rev(spec@presample[, 3]), rep(NA, times = n))
    m = length(spec@presample[, 1])
    names(z) = names(h) = names(y) = NULL
    mu = model$mu
    ar = model$ar
    ma = model$ma
    omega = model$omega
    alpha = model$alpha
    gamma = model$gamma
    beta = model$beta
    deltainv = 1/delta
    order.ar = length(ar)
    order.ma = length(ma)
    order.alpha = length(alpha)
    order.beta = length(beta)
    eps = h^deltainv * z
    for (i in (m + 1):(n + m)) {
        h[i] = omega + sum(alpha * (abs(eps[i - (1:order.alpha)]) - 
            gamma * (eps[i - (1:order.alpha)]))^delta) + sum(beta * 
            h[i - (1:order.beta)])
        eps[i] = h[i]^deltainv * z[i]
        y[i] = sum(ar * y[i - (1:order.ar)]) + sum(ma * eps[i - 
            (1:order.ma)]) + eps[i]
    }
    y = y + mu
    data = cbind(z = z[(m + 1):(n + m)], sigma = h[(m + 1):(n + 
        m)]^deltainv, y = y[(m + 1):(n + m)])
    rownames(data) = as.character(1:n)
    if (n.start > 0) 
        data = data[-(1:n.start), ]
    from <- timeDate(format(Sys.time(), format = "\%Y-\%m-\%d")) - 
        NROW(data) * 24 * 3600
    charvec <- timeSequence(from = from, length.out = NROW(data))
    ans <- timeSeries(data = data[, c(3, 2, 1)], charvec = charvec)
    colnames(ans) <- c("Series", "Volatility", "Innovations")
    attr(ans, "control") <- list(garchSpec = spec)
    ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
