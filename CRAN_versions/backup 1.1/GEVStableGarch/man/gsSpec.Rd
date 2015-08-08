\name{gsSpec}


\alias{gsSpec}


\title{Specification of ARMA-GARCH/APARCH models with GEV or stable distributions}


\description{

Specifies an ARMA-GARCH or ARMA-APARCH model with innovations following GEV, stable.

}
\usage{
gsSpec(model = list(), presample = NULL, 
    cond.dist = c("stableS0", "stableS1", "stableS2", "gev", 
    "gat", "norm", "std", "sstd", "skstd", "ged"), 
    rseed = NULL)
}


\arguments{

    \item{model}{
        a list of ARMA-GARCH/APARCH model parameters: \cr
        \code{omega} - the coefficient of the variance equation; \cr
        \code{alpha} - the value or vector of autoregressive coefficients; \cr
        \code{beta} - the value or vector of variance coefficients; 
        \cr
        The values for the ARMA part are: \cr
        \code{mu} - the mean value; \cr
        \code{ar} - the autoregressive coefficients; \cr
        \code{ma} - the moving average coefficients. 
        \cr
        The parameters for the conditional distributions are:\cr
        \code{skew} - the skewness parameter; \cr
        \code{shape} - the shape parameter.\cr
        }   

    \item{presample}{
        presample - a numeric "matrix" with 3 columns and
        with max(m,n,p,q) rows. 
        The first culumn are the innovations, the second
        the conditional variances, and the last the time series.
        When presample is missing, it is constructed
        as [z,h,y] where z ~ Normal(0,1), h = "uev" recursion 
        initialization described in Wuertz et al. (2009) and
        y = mu. Note that the conditional variance column
        can contain only strictly positive numbers and the function
        \code{gsSpec} check for invalid values.
        }
        
        \item{cond.dist}{ a character string naming conditional distribution of innovations. The package was created to accept the following distributions: \code{"stableS0"}, \code{"stableS1"}, \code{"stableS2"}, \code{"gev"} and \code{"gat"}. Other common distributions are also possible such as the \code{"norm"}, \code{"std"}, \code{"sstd"}, \code{"skstd"}, \code{"ged"}, but they are not the main feature of this work.
        }
        
    \item{rseed}{
        the seed for the intitialization of the random number generator for the innovations.
        }
}

\details{

This functions uses the interface of the \code{garchSpec} routine from package \pkg{fGarch} to simulate random values of the ARMA-GARCH/APARCH model with conditional GEV or stable distribution. 

}

\value{

    The returned value is an object of class \code{"GEVSTABLEGARCHSPEC"}.
    
}

\references{

    Wuertz, D., Chalabi, Y., with contribution from Miklovic, M., Boudt, C., Chausse, P., and others (2013).
    \emph{fGarch: Rmetrics - Autoregressive Conditional Heteroskedastic Modelling, R package version 3010.82},   
    \url{http://CRAN.R-project.org/package=fGarch}.

    Wuertz, D., Chalabi, Y., Luksan, L. (2009).
    \emph{Parameter Estimation of ARMA Models with GARCH/ APARCH Errors: An R and SPlus SoftwareImplementation}.   
    Journal of Statistical Software, forthcoming,
    \href{http://www-stat.wharton.upenn.edu/~steele/Courses/956/RResources/GarchAndR/WurtzEtAlGarch.pdf}{http://www-stat.wharton.upenn.edu/~steele/...WurtzEtAlGarch.pdf}.

}

\author{

Thiago do Rego Sousa.

}

\examples{

# stable-GARCH from Curto et al. (2009) for the DJIA dataset
spec.stable = gsSpec(model = list(mu = 0.0596, omega = 0.0061, 
alpha = 0.0497, beta = 0.9325, skew = -0.9516, shape = 1.9252), 
cond.dist = "stableS1")
sim.stable = gsSim(spec = spec.stable, n = 1000)
 
# GEV-GARCH model from Zhao et al. (2011)
spec.gev = gsSpec(model = list(mu = 0.21, a = 0.32, omega = 0.01,
alpha = 0.45, beta = 0.08, shape = 0.08), cond.dist = "gev")
sim.gev = gsSim(spec = spec.gev, n = 1000)

}

\keyword{ models }
