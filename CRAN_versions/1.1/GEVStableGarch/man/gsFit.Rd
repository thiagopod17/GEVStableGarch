\name{gsFit}
\alias{gsFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsFit(formula = ~garch(1, 1), data, cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "GAt", "norm", "std", "sstd", "skstd", "ged"), include.mean = TRUE, algorithm = c("sqp", "sqp.restriction", "nlminb", "nlminb+nm"), control = NULL, tolerance = NULL, title = NULL, description = NULL, DEBUG = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
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
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{tolerance}{
%%     ~~Describe \code{tolerance} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{description}{
%%     ~~Describe \code{description} here~~
}
  \item{DEBUG}{
%%     ~~Describe \code{DEBUG} here~~
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
function (formula = ~garch(1, 1), data, cond.dist = c("stableS0", 
    "stableS1", "stableS2", "gev", "GAt", "norm", "std", "sstd", 
    "skstd", "ged"), include.mean = TRUE, algorithm = c("sqp", 
    "sqp.restriction", "nlminb", "nlminb+nm"), control = NULL, 
    tolerance = NULL, title = NULL, description = NULL, DEBUG = FALSE) 
{
    cond.dist = match.arg(cond.dist)
    algorithm = match.arg(algorithm)
    if (!is.numeric(data) || any(is.na(data)) || any(is.null(data)) || 
        any(!is.finite(data))) 
        stop("Invalid 'data' input. It may be contain NA, NULL or Inf.")
    CALL = match.call()
    if (is.null(tolerance)) 
        tolerance = list(TOLG = 1e-08, TOLSTABLE = 0.01, TOLSTATIONARITY = 0.001)
    if (tolerance$TOLSTATIONARITY > 0.05) 
        stop("TOLSTATIONARITY can not be > 0.05")
    if (is.null(title)) 
        title = "ARMA-GARCH modelling"
    formula.input <- formula
    formula <- .getFormula(formula)
    m <- formula$formula.order[1]
    n <- formula$formula.order[2]
    p <- formula$formula.order[3]
    q <- formula$formula.order[4]
    APARCH <- formula$isAPARCH
    formula.mean <- ""
    formula.var <- ""
    if (m > 0 || n > 0) 
        formula.mean <- formula$formula.mean
    else formula.mean <- ""
    if (p > 0) 
        formula.var <- formula$formula.var
    else formula.var <- ""
    if (algorithm == "sqp.restriction") {
        if (formula.var == "") 
            stop("sqp.restriction should only be used with GARCH/APARCH models.")
        if (any(cond.dist == c("stableS0", "stableS2"))) 
            stop("sqp.restriction algorithm can only be used with \n stable distribution in S1 parametrization, i.e., cond.dist = 'stableS1'.")
    }
    if (formula.var == "") 
        stop("Pure ARMA model not allowed")
    printRes = TRUE
    ARMAonly <- FALSE
    AR <- FALSE
    MA <- FALSE
    GARCH <- FALSE
    if (m == 0) 
        AR <- TRUE
    if (n == 0) 
        MA <- TRUE
    if (q == 0) 
        GARCH <- TRUE
    optim.finished <- FALSE
    if (AR == TRUE) 
        m <- 1
    if (MA == TRUE) 
        n <- 1
    if ((p == 0) && (q == 0)) 
        ARMAonly = TRUE
    if (GARCH == TRUE && !ARMAonly) 
        q <- 1
    data <- data
    N <- length(data)
    out <- NULL
    messages <- NULL
    out$order <- c(m, n, p, q, include.mean, APARCH)
    TMPvector <- c(if (m != 0 || n != 0) c("arma(", m, ",", n, 
        ")-"), if (APARCH == FALSE) c("garch(", p, ",", q, ")"), 
        if (APARCH == TRUE) c("aparch(", p, ",", q, ")"))
    TMPorder <- paste(TMPvector, collapse = "")
    TMPvectorintercept <- c("include.mean:", if (include.mean == 
        TRUE) "TRUE", if (include.mean == FALSE) "FALSE")
    TMPintercept <- paste(TMPvectorintercept, collapse = "")
    out$model <- paste(TMPorder, "##", TMPintercept, collapse = "")
    out$cond.dist <- cond.dist
    out$data <- data
    optim.finished <- FALSE
    if (cond.dist == "GAt") 
        lengthShape = 2
    else lengthShape = 1
    if (DEBUG) 
        print(c("lengthShape", lengthShape))
    arCheck <- function(ar) {
        p <- max(which(c(1, -ar) != 0)) - 1
        if (!p) 
            return(TRUE)
        all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
    }
    armaLLH <- function(parm) {
        if (DEBUG) 
            print(parm)
        if (sum(is.nan(parm)) != 0) {
            return(1e+99)
        }
        mu <- parm[1]
        a <- parm[(1 + 1):(2 + m - 1)]
        b <- parm[(1 + m + 1):(2 + m + n - 1)]
        skew <- parm[1 + m + n + 1]
        shape <- parm[(2 + m + n + 1):(2 + m + n + lengthShape)]
        sigma <- parm[2 + lengthShape + m + n + 1]
        if (AR == TRUE) 
            a <- 0
        if (MA == TRUE) 
            b <- 0
        if (include.mean == FALSE) 
            mu <- 0
        if (any(cond.dist == c("stableS0", "stableS1", "stableS2"))) {
            if (shape - delta < tolerance$TOLSTABLE || abs(shape) < 
                tolerance$TOLSTABLE || !(abs(skew) < 1) || !((shape - 
                2) < 0)) {
                return(1e+99)
            }
        }
        if (!arCheck(a)) 
            return(1e+99)
        if (sigma <= 0) 
            return(1e+99)
        z <- .filterArma(data = data, size.data = N, m = m, n = n, 
            mu = mu, a = a, b = b)
        if (DEBUG) {
            print(c("length(z)", length(z)))
            print(c("N", N))
        }
        if (optim.finished) {
            out$residuals <<- as.numeric(z)
            out$sigma.t <<- sigma
            out$h.t <<- sigma
        }
        llh.dens <- .armaDist(z = z, sigma = sigma, shape = shape, 
            skew = skew, cond.dist = cond.dist)
        llh <- llh.dens
        if (is.nan(llh) || is.infinite(llh) || is.na(llh)) {
            llh <- 1e+99
        }
        llh
    }
    garchLLH = function(parm) {
        if (sum(is.nan(parm)) != 0) {
            return(1e+99)
        }
        mu <- parm[1]
        a <- parm[(1 + 1):(2 + m - 1)]
        b <- parm[(1 + m + 1):(2 + m + n - 1)]
        omega <- parm[1 + m + n + 1]
        alpha <- parm[(2 + m + n + 1):(3 + m + n + p - 1)]
        gm <- parm[(2 + m + n + p + 1):(3 + m + n + p + p - 1)]
        beta <- parm[(2 + m + n + 2 * p + 1):(3 + m + n + 2 * 
            p + q - 1)]
        delta <- parm[2 + m + n + 2 * p + q + 1]
        skew <- parm[3 + m + n + 2 * p + q + 1]
        shape <- parm[(4 + m + n + 2 * p + q + 1):(4 + m + n + 
            2 * p + q + lengthShape)]
        if (!APARCH) {
            gm = rep(0, p)
            delta = 2
            if (any(cond.dist == c("stableS0", "stableS1", "stableS2"))) 
                delta = 1
        }
        if (AR == TRUE) 
            a <- 0
        if (MA == TRUE) 
            b <- 0
        if (GARCH == TRUE) 
            beta <- 0
        if (include.mean == FALSE) 
            mu <- 0
        cond.general <- FALSE
        cond.normal <- FALSE
        cond.student <- FALSE
        cond.gev <- FALSE
        cond.stable <- FALSE
        parset <- c(omega, alpha, if (!GARCH) beta, delta)
        cond.general <- any(parset < tolerance$TOLG)
        if (cond.general || cond.student || cond.gev || cond.stable) {
            return(1e+99)
        }
        if (!arCheck(a)) 
            return(1e+99)
        if (AR == TRUE && MA == TRUE) {
            filteredSeries <- .filterAparch(data = data, p = p, 
                q = q, mu = mu, omega = omega, alpha = alpha, 
                beta = beta, gamma = gm, delta = delta)
            z <- filteredSeries[, 1]
            hh <- filteredSeries[, 2]
        }
        if (AR == FALSE || MA == FALSE) {
            filteredArma <- .filterArma(data = data, size.data = N, 
                m = m, n = n, mu = mu, a = a, b = b)
            filteredSeries <- .filterAparch(data = filteredArma, 
                p = p, q = q, mu = 0, omega = omega, alpha = alpha, 
                beta = beta, gamma = gm, delta = delta)
            z <- filteredSeries[, 1]
            hh <- filteredSeries[, 2]
        }
        if (optim.finished) {
            out$residuals <<- as.numeric(z)
            out$sigma.t <<- as.numeric(hh)
            out$h.t <<- as.numeric(hh^delta)
        }
        llh.dens <- .armaGarchDist(z = z, hh = hh, shape = shape, 
            skew = skew, cond.dist = cond.dist)
        llh <- llh.dens
        if (is.nan(llh) || is.infinite(llh) || is.na(llh)) {
            llh <- 1e+99
        }
        llh
    }
    start <- .getStart(data = data, m = m, n = n, p = p, q = q, 
        AR = AR, MA = MA, cond.dist = cond.dist, TOLG = tolerance$TOLG, 
        TOLSTABLE = tolerance$TOLSTABLE)
    if (DEBUG) {
        print("start")
        print(start)
    }
    garch.stationarity <- function(parm) {
        omega <- parm[1 + m + n + 1]
        alpha <- parm[(2 + m + n + 1):(3 + m + n + p - 1)]
        gm <- parm[(2 + m + n + p + 1):(3 + m + n + p + p - 1)]
        beta <- parm[(2 + m + n + 2 * p + 1):(3 + m + n + 2 * 
            p + q - 1)]
        delta <- parm[2 + m + n + 2 * p + q + 1]
        skew <- parm[3 + m + n + 2 * p + q + 1]
        shape <- parm[(4 + m + n + 2 * p + q + 1):(4 + m + n + 
            2 * p + q + lengthShape)]
        model = list(omega = omega, alpha = alpha, gm = gm, beta = beta, 
            delta = delta, skew = skew, shape = shape)
        .stationarityAparch(model = model, formula = formula, 
            cond.dist = cond.dist)
    }
    if (ARMAonly) 
        modelLLH <- armaLLH
    else modelLLH <- garchLLH
    if (algorithm == "nlminb") {
        fit1 <- nlminb(start[1, ], objective = modelLLH, lower = start[2, 
            ], upper = start[3, ], control = control)
        out$llh <- fit1$objective
        out$par <- fit1$par
        out$hessian <- fit1$hessian
        out$convergence <- fit1$convergence
    }
    if (algorithm == "nlminb+nm") {
        fit1.partial <- nlminb(start[1, ], objective = modelLLH, 
            lower = start[2, ], upper = start[3, ], control = control)
        fit1 <- optim(par = fit1.partial$par, fn = modelLLH, 
            method = "Nelder-Mead", hessian = TRUE)
        out$llh <- fit1$value
        out$par <- fit1$par
        out$hessian <- fit1$hessian
        out$convergence <- fit1$convergence
    }
    if (algorithm == "sqp") {
        fit1 <- solnp(pars = start[1, ], fun = modelLLH, LB = start[2, 
            ], UB = start[3, ], control = control)
        out$llh <- fit1$values[length(fit1$values)]
        out$par <- fit1$pars
        out$hessian <- fit1$hessian
        out$convergence <- fit1$convergence
    }
    if (algorithm == "sqp.restriction") {
        fit1 <- solnp(pars = start[1, ], fun = modelLLH, ineqfun = garch.stationarity, 
            ineqLB = 0, ineqUB = 1 - tolerance$TOLSTATIONARITY, 
            LB = start[2, ], UB = start[3, ], control = control)
        out$llh <- fit1$values[length(fit1$values)]
        out$par <- fit1$pars
        out$hessian <- fit1$hessian
        sizeHessian = length(out$hessian[1, ]) - 1
        out$hessian = out$hessian[1:sizeHessian, 1:sizeHessian]
        out$convergence <- fit1$convergence
    }
    if ((algorithm == "sqp") && !is.numeric(try(sqrt(diag(solve(out$hessian))), 
        silent = TRUE))) {
        fit1.partial <- optim(par = fit1$par, fn = modelLLH, 
            method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 1))
        if (sum(abs(fit1.partial$par - fit1$par)) == 0) {
            out$hessian = fit1.partial$hessian
            if (DEBUG) 
                print(abs(fit1.partial$par - fit1$par))
        }
    }
    if ((algorithm == "sqp.restriction")) {
        fit1.partial <- optim(par = fit1$par, fn = modelLLH, 
            method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 1))
        if (sum(abs(fit1.partial$par - fit1$par)) == 0) {
            out$hessian = fit1.partial$hessian
            if (DEBUG) 
                print(abs(fit1.partial$par - fit1$par))
        }
    }
    if (DEBUG) 
        print(fit1)
    if (out$llh == 1e+99) 
        out$convergence = 1
    optim.finished = TRUE
    modelLLH(out$par)
    if (!ARMAonly) {
        outindex <- c(if (include.mean) 1, if (AR == FALSE) (1 + 
            1):(2 + m - 1), if (MA == FALSE) (1 + m + 1):(2 + 
            m + n - 1), if (!ARMAonly) (1 + m + n + 1), if (!ARMAonly) (2 + 
            m + n + 1):(3 + m + n + p - 1), if (APARCH) (2 + 
            m + n + p + 1):(3 + m + n + p + p - 1), if (!GARCH) (2 + 
            m + n + 2 * p + 1):(3 + m + n + 2 * p + q - 1), if (APARCH) (2 + 
            m + n + 2 * p + q + 1), if (any(c("sstd", "skstd", 
            "stableS0", "stableS1", "stableS2", "GAt") == cond.dist)) (3 + 
            m + n + 2 * p + q + 1), if (any(c("std", "gev", "stableS0", 
            "stableS1", "stableS2", "sstd", "skstd", "ged", "GAt") == 
            cond.dist)) (4 + m + n + 2 * p + q + 1):(4 + m + 
            n + 2 * p + q + lengthShape))
    }
    else {
        outindex <- c(if (include.mean) 1, if (AR == FALSE) (1 + 
            1):(2 + m - 1), if (MA == FALSE) (1 + m + 1):(2 + 
            m + n - 1), if (!ARMAonly) (1 + m + n + 1), if (!ARMAonly) (2 + 
            m + n + 1):(3 + m + n + p - 1), if (APARCH) (2 + 
            m + n + p + 1):(3 + m + n + p + p - 1), if (!GARCH) (2 + 
            m + n + 2 * p + 1):(3 + m + n + 2 * p + q - 1), if (APARCH) (2 + 
            m + n + 2 * p + q + 1), if (any(c("sstd", "skstd", 
            "stableS0", "stableS1", "stableS2", "GAt") == cond.dist)) (1 + 
            m + n + 2 * p + q + 1), if (any(c("std", "gev", "stableS0", 
            "stableS1", "stableS2", "sstd", "skstd", "ged", "GAt") == 
            cond.dist)) (2 + m + n + 2 * p + q + 1):(2 + m + 
            n + 2 * p + q + lengthShape), length(out$par))
    }
    outnames <- c(if (include.mean) "mu", if (AR == FALSE) paste("ar", 
        1:m, sep = ""), if (MA == FALSE) paste("ma", 1:n, sep = ""), 
        if (!ARMAonly) "omega", if (!ARMAonly) paste("alpha", 
            1:p, sep = ""), if (APARCH) paste("gamma", 1:p, sep = ""), 
        if (!GARCH) paste("beta", 1:q, sep = ""), if (APARCH) "delta", 
        if (any(c("sstd", "stableS0", "stableS1", "stableS2", 
            "GAt", "skstd") == cond.dist)) "skew", if (any(c("std", 
            "gev", "stableS0", "stableS1", "stableS2", "sstd", 
            "ged", "GAt", "skstd") == cond.dist)) paste("shape", 
            1:lengthShape, sep = ""), if (ARMAonly) "sigma")
    if (DEBUG) {
        print(c("out", out))
        print(c("outindex", outindex))
    }
    out$par <- out$par[outindex]
    names(out$par) <- outnames
    out$hessian <- out$hessian[outindex, outindex]
    nParam <- length(out$par)
    out$aic = 2 * out$llh + 2 * nParam
    out$aicc = 2 * out$llh + 2 * (nParam + 1) * N/(N - nParam - 
        2)
    out$bic = 2 * out$llh + nParam * log(N)
    out$ics = c(out$aic, out$bic, out$aicc)
    names(out$ics) <- c("AIC", "BIC", "AICc")
    if (printRes) {
        solveHessianFailed = FALSE
        out$se.coef <- 0
        out$se.coef <- try(sqrt(diag(solve(out$hessian))), silent = TRUE)
        if (!is.numeric(out$se.coef)) {
            solveHessianFailed = TRUE
            messages$solving.hessian.matrix = "Error inverting Hessian Matrix"
            out$matcoef <- cbind(out$par, rep(NA, length(out$par)), 
                rep(NA, length(out$par)), rep(NA, length(out$par)))
            dimnames(out$matcoef) = dimnames(out$matcoef) = list(names(out$par), 
                c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
        }
        else {
            out$tval <- try(out$par/out$se.coef, silent = TRUE)
            out$matcoef = cbind(out$par, if (is.numeric(out$se.coef)) 
                out$se.coef, if (is.numeric(out$tval)) 
                out$tval, if (is.numeric(out$tval)) 
                2 * (1 - pnorm(abs(out$tval))))
            dimnames(out$matcoef) = list(names(out$tval), c(" Estimate", 
                " Std. Error", " t value", "Pr(>|t|)"))
        }
        cat("\nFinal Estimate of the Negative LLH:\n")
        cat("-LLH:", out$llh)
        if (out$convergence == 0) 
            messages$optimization.algorithm = "Algorithm achieved convergence"
        else messages$optimization.algorithm = "Algorithm did not achieved convergence"
        cat("\nCoefficient(s):\n")
        printCoefmat(round(out$matcoef, digits = 6), digits = 6, 
            signif.stars = TRUE)
    }
    out$order <- c(formula$formula.order[1], formula$formula.order[2], 
        formula$formula.order[3], formula$formula.order[4])
    names(out$order) <- c("m", "n", "p", "q")
    fit <- list(par = out$par, llh = out$llh, hessian = out$hessian, 
        ics = out$ics, order = out$order, cond.dist = cond.dist, 
        se.coef = out$se.coef, tval = out$tval, matcoef = out$matcoef)
    new("fGEVSTABLEGARCH", call = as.call(match.call()), formula = formula.input, 
        method = "Max Log-Likelihood Estimation", convergence = out$convergence, 
        messages = messages, data = data, fit = fit, residuals = out$residuals, 
        h.t = out$h.t, sigma.t = as.vector(out$sigma.t), title = as.character(title), 
        description = as.character(description))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
