\name{gsMomentAparch}
\alias{gsMomentAparch}


\title{

Computation of moments for several conditional distribution

}


\description{

Computation of the moments expression E( |Z| - gamma * Z) ^ delta

}


\usage{

gsMomentAparch(cond.dist = c("stableS1", "gev", "GAt", "norm", 
"std", "sstd", "skstd", "ged"), shape = 1.5, skew = 0, delta = 1, gm = 0)

}


\arguments{

        \item{cond.dist}{ a character string naming conditional distribution of innovations. The package was created to accept the following distributions: \code{"stableS0"}, \code{"stableS1"}, \code{"stableS2"}, \code{"gev"} and \code{"GAt"}. Other common distributions are also possible such as the \code{"norm"}, \code{"std"}, \code{"sstd"}, \code{"skstd"}, \code{"ged"}, but they are not the main feature of this work.
        }
        
  \item{shape, skew}{
        The shape and skew parameter for the conditional distribution.
}

  \item{delta, gm}{
        The delta and gm parameteres of the ARMA-GARCH/APARCH model
}

}


\details{
EXPLICAR COMO EH FEITO ESSES CALCULOS. COLOCAR EXPRESSOES DO ARTIGO.
}


\value{

Returns the following expression for several conditional distributions
E( |Z| - gamma * Z) ^ delta

}


\references{

NOSSO ARTIGO JSS

Diongue AK (2008). "An investigation of Stable-Paretian Asymmetric Power GARCH Model." Journal des Sciences, 8(4), 15-26.

Mittnik S, Paolella MS (2000). "Conditional Density and Value-At-Risk Prediction of Asian Currency Exchange Rates." Journal of Forecasting, 19(4), 313-333.

Mittnik S, Paolella M, Rachev S (2002). "Stationarity of Stable Power-GARCH Processes." Journal of Econometrics, 106, 97-107.

}


\author{

Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano 
and Silvia Regina Costa Lopes

}


\examples{

# EXEMPLOS GSMOMENTAPARCH

}


\keyword{ aparch }
\keyword{ moments }