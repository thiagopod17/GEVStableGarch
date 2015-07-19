\name{gsSpec}
\alias{gsSpec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsSpec(model = list(), presample = NULL, cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "GAt", "norm", "std", "sstd", "skstd", "ged"), rseed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{presample}{
%%     ~~Describe \code{presample} here~~
}
  \item{cond.dist}{
%%     ~~Describe \code{cond.dist} here~~
}
  \item{rseed}{
%%     ~~Describe \code{rseed} here~~
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
function (model = list(), presample = NULL, cond.dist = c("stableS0", 
    "stableS1", "stableS2", "gev", "GAt", "norm", "std", "sstd", 
    "skstd", "ged"), rseed = NULL) 
{
    skew = list(stableS0 = 0, stableS1 = 0, stableS2 = 0, gev = NULL, 
        GAt = 1, norm = NULL, std = NULL, sstd = 0.9, skstd = 1, 
        ged = NULL)
    shape = list(stableS0 = 1.7, stableS1 = 1.7, stableS2 = 1.7, 
        gev = 0.3, GAt = c(3, 1), norm = NULL, std = 4, sstd = 4, 
        skstd = 4, ged = 2)
    cond.dist = match.arg(cond.dist)
    initialDelta = NULL
    if (!is.null(model$alpha) && is.null(model$delta)) {
        if (any(cond.dist == c("stableS0", "stableS1", "stableS2"))) 
            initialDelta = 1
        else initialDelta = 2
    }
    control = list(omega = 1, alpha = NULL, gamma = NULL, beta = NULL, 
        mu = NULL, ar = NULL, ma = NULL, delta = initialDelta, 
        skew = skew[[cond.dist]], shape = shape[[cond.dist]])
    control[names(model)] <- model
    model <- control
    order.ar = length(model$ar)
    order.ma = length(model$ma)
    order.alpha = length(model$alpha)
    order.beta = length(model$beta)
    order.max = max(order.ar, order.ma, order.alpha, order.beta)
    if ((order.alpha == 0 && order.beta != 0)) 
        stop("In a Garch(p,q)/Aparch(p,q) model we must have p > 0")
    if (!is.null(model$delta)) {
        if (length(model$delta) > 1) 
            stop("The parameter 'delta' must be a single number")
        if (!(model$delta > 0)) 
            stop("The parameter 'delta' must be > 0.")
    }
    if ((length(model$gamma) != 0) && (length(model$alpha) != 
        0)) {
        if (length(model$alpha) != length(model$gamma)) 
            stop("'alpha' and 'gamma' must have the same size for APARCH models")
    }
    if (!is.null(model$gamma)) {
        if (sum(!(abs(model$gamma) < 1)) > 0) 
            stop("The parameter 'gamma' must be in the range -1 < gamma < 1")
    }
    if (sum(model$alpha < 0) > 0) 
        stop("The parameter 'alpha' must be greater than or equal 0")
    if (sum(model$beta < 0) > 0) 
        stop("The parameter 'beta' must be greater than or equal 0")
    if (!(model$omega > 0)) 
        stop("The parameter 'omega' must be > 0")
    if (is.null(model$ar) && is.null(model$ma) && is.null(model$alpha)) 
        stop("The model parameters were not specified correctly.")
    if (!is.null(model$delta) && is.null(model$alpha)) 
        stop("The parameter delta should be only specified for a GARCH or APARCH model.")
    if (is.matrix(presample) && !is.null(presample)) {
        if (dim(presample)[2] != 3 || dim(presample)[1] < order.max) 
            stop(cat("The presample object should be a matrix with three columns formated as: \n [Innovations, Conditional Variance, Time Series] with dimensions \n l x 3, where l = max(m,n,p,q). "))
        if (dim(presample)[1] != order.max) {
            warning(cat("The number of columns of the Presample matrix is \n bigger than l = max(m,n,p,q). The simulated series \n will use only the first 'l' columns"))
            presample = as.matrix(presample[1:order.max, ])
            if (order.max == 1) 
                presample = t(presample)
        }
    }
    formula.mean = ""
    if (order.ar == 0 && order.ma == 0) {
        formula.mean = ""
    }
    else {
        formula.mean = paste("arma(", as.character(order.ar), 
            ", ", as.character(order.ma), ")", sep = "")
    }
    formula.var = ""
    if (order.alpha > 0) 
        formula.var = "garch"
    if (!is.null(model$alpha)) {
        if (!is.null(model$gamma)) {
            if (sum(model$gamma == 0) != length(model$gamma)) 
                formula.var = "aparch"
        }
        else {
            if (model$delta != 2 && !any(cond.dist == c("stableS0", 
                "stableS1", "stableS2"))) 
                formula.var = "aparch"
            if (model$delta != 1 && any(cond.dist == c("stableS0", 
                "stableS1", "stableS2"))) 
                formula.var = "aparch"
        }
    }
    if (order.alpha == 0 && order.beta == 0) {
        formula.var = formula.var
    }
    if (order.alpha > 0) {
        formula.var = paste(formula.var, "(", as.character(order.alpha), 
            ", ", as.character(order.beta), ")", sep = "")
    }
    if (formula.mean == "") {
        formula = as.formula(paste("~", formula.var))
    }
    if (formula.var == "") {
        formula = as.formula(paste("~", formula.mean))
    }
    if ((formula.mean != "") && (formula.var != "")) {
        formula = as.formula(paste("~", formula.mean, "+", formula.var))
    }
    if (formula.var == "") 
        stop("Pure ARMA model not allowed")
    if (is.null(model$mu)) 
        model$mu = 0
    if (is.null(model$ar)) 
        model$ar = 0
    if (is.null(model$ma)) 
        model$ma = 0
    if (is.null(model$gamma)) 
        model$gamma = rep(0, times = order.alpha)
    if (is.null(rseed)) {
        rseed = 0
    }
    else {
        set.seed(rseed)
    }
    persistency = 1 - sum(model$alpha) - sum(model$beta)
    if (persistency * (1 - persistency) < 0) 
        persistency = 0.1
    if (!is.null(presample) && is.matrix(presample)) {
        z = presample[, 1]
        h = presample[, 2]
        y = presample[, 3]
        if (sum(!(h > 0))) 
            stop("Conditional Variance column can have only strictly positive numbers")
    }
    else {
        z = rnorm(n = order.max)
        h = rep(model$omega * (1 + persistency * (1 - persistency)), 
            times = order.max)
        y = rep(model$mu, times = order.max)
    }
    presample = cbind(z, h, y)
    new("fGEVSTABLEGARCHSPEC", call = match.call(), formula = formula, 
        model = list(omega = model$omega, alpha = model$alpha, 
            gamma = model$gamma, beta = model$beta, mu = model$mu, 
            ar = model$ar, ma = model$ma, delta = model$delta, 
            skew = model$skew, shape = model$shape), presample = as.matrix(presample), 
        distribution = as.character(cond.dist), rseed = as.numeric(rseed))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
