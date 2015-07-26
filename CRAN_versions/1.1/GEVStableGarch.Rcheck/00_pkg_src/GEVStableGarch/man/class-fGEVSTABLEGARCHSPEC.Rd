\name{fGEVSTABLEGARCHSPEC-class}


\Rdversion{1.1}


\docType{class}


\alias{fGEVSTABLEGARCHSPEC-class}


\alias{show,fGEVSTABLEGARCHSPEC-update}


\title{Class \code{"fGEVSTABLEGARCHSPEC"}}


\description{

This class represents an ARMA-GARCH/APARCH model with innovations following an i.i.d
stable, GEV distributions. Other distributions are also allowed.

}


\section{Objects from the Class}{

Objects can be created by calls of the function \code{gsSpec}. This object 
specifies an univariate ARMA-GARCH/APARCH model with GEV or stable conditional
distribution. Note that we only allow the user
to create an object that do not violate the model specifications. For more details about
the definition of those models see:
see Mittnik et al. (2002) (stable innovations), Zhao et al. (2011) (GEV innovations) and Wuertz et al. (2008) (finite variance innovations).

}


\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"call"}: 
        the call of the \code{gsSpec} function.
        } 
    \item{\code{formula}:}{Object of class \code{"formula"}: 
        a list with two formula entries for the mean and variance 
        equation of the combined ARMA-GARCH/APARCH model.
        }
    \item{\code{model}:}{Object of class \code{"list"}: 
        a list with the model parameters specifying the Autorregressive,
        the GARCH/APARCH and the innovations parameters.
        }
    \item{\code{presample}:}{Object of class \code{"matrix"}: 
        a numeric matrix with presample values.
        }
    \item{\code{distribution}:}{Object of class \code{"character"}: 
        a character string with the name of the conditional distribution.
        }
    \item{\code{rseed}:}{Object of class \code{"numeric"}: 
        an integer with the random number generator seed.
        }
  }
}


\section{Methods}{
  \describe{
  
    \item{show}{\code{signature(object = "fGEVSTABLEGARCHSPEC")}: 
    prints an object of class 'fGEVSTABLEGARCH' 
    }
    
	 }
}


\references{

Mittnik S, Paolella M, Rachev S (2002). "Stationarity of Stable Power-GARCH Processes." Journal of Econometrics, 106, 97-107.

Wuertz D, Chalabi Y, Luksan L (2009). "Parameter Estimation of ARMA Models with GARCH/APARCH Errors: An R and SPlus Software Implementation." Journal of Statistical Software, forthcoming. URL http://www-stat.wharton.upenn.edu/~steele/Courses/956/RResources/GarchAndR/WurtzEtAlGarch.pdf.

Zhao X, Scarrott C, Oxley L, Reale M (2011). "GARCH Dependence in Extreme Value Models with Bayesian Inference." Mathematics and Computers in Simulation, 81(7), 1430-1440.

}


\author{

Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano 
and Silvia Regina Costa Lopes.

}


\keyword{classes}
