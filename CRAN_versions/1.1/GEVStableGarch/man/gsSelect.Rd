\name{gsSelect}
\alias{gsSelect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsSelect(data, order.max = c(1, 1, 1, 1), selection.criteria = c("AIC", "AICc", "BIC"), is.aparch = FALSE, cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "GAt", "norm", "std", "sstd", "skstd", "ged"), include.mean = TRUE, algorithm = c("sqp", "sqp.restriction", "nlminb", "nlminb+nm"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{order.max}{
%%     ~~Describe \code{order.max} here~~
}
  \item{selection.criteria}{
%%     ~~Describe \code{selection.criteria} here~~
}
  \item{is.aparch}{
%%     ~~Describe \code{is.aparch} here~~
}
  \item{cond.dist}{
%%     ~~Describe \code{cond.dist} here~~
}
  \item{include.mean}{
%%     ~~Describe \code{include.mean} here~~
}
  \item{algorithm}{
%%     ~~Describe \code{algorithm} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (data, order.max = c(1, 1, 1, 1), selection.criteria = c("AIC", 
    "AICc", "BIC"), is.aparch = FALSE, cond.dist = c("stableS0", 
    "stableS1", "stableS2", "gev", "GAt", "norm", "std", "sstd", 
    "skstd", "ged"), include.mean = TRUE, algorithm = c("sqp", 
    "sqp.restriction", "nlminb", "nlminb+nm"), ...) 
{
    cond.dist = match.arg(cond.dist)
    algorithm = match.arg(algorithm)
    selection.criteria = match.arg(selection.criteria)
    if (!is.numeric(data) || !is.vector(data)) 
        stop("data set must be a numerical one dimensional vector")
    goodness.of.fit.min <- 1e+99
    goodness.of.fit.current <- 1e+99
    fit.min <- list()
    fit <- list()
    order.list <- .getOrder(order.max = order.max)
    order.list.size <- length(order.list[, 1])
    for (i in 1:order.list.size) {
        m = order.list[i, 1]
        n = order.list[i, 2]
        p = order.list[i, 3]
        q = order.list[i, 4]
        if (is.aparch == TRUE) {
            formula = as.formula(paste("~ arma(", m, ", ", n, 
                ") + aparch(", p, ", ", q, ")", sep = "", collapse = NULL))
        }
        else {
            formula = as.formula(paste("~ arma(", m, ", ", n, 
                ") + garch(", p, ", ", q, ")", sep = "", collapse = NULL))
        }
        cat("\n------------------------------------------------------------------------------------------\n")
        cat(paste("Model: ", paste(formula)[2], " with '", cond.dist, 
            "' conditional distribution", sep = ""))
        cat("\n------------------------------------------------------------------------------------------\n")
        fit = gsFit(data = data, formula = formula, cond.dist = cond.dist, 
            algorithm = algorithm, ...)
        goodness.of.fit.current = fit@fit$ics[selection.criteria]
        if (goodness.of.fit.current < goodness.of.fit.min) {
            fit.min <- fit
            goodness.of.fit.min <- fit@fit$ics[selection.criteria]
        }
    }
    cat("\n------------------------------------------------------------------------------------------\n")
    cat(paste("Best Model: ", paste(fit.min@formula)[2]))
    cat("\n------------------------------------------------------------------------------------------\n")
    return(fit.min)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
