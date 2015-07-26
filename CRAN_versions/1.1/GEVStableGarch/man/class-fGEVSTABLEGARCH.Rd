\name{fGEVSTABLEGARCH-class}


\Rdversion{1.1}


\docType{class}


\alias{fGEVSTABLEGARCH-class}


\alias{show,fGEVSTABLEGARCH-update}


\title{Class \code{"fGEVSTABLEGARCH"}}


\description{

The class fGEVSTABLEGARCH represents an ARMA-GARCH/APARCH 
model estimated wih function \code{gsFit}

}


\section{Objects from the Class}{

Objects can be created by calls of the form \code{gsFit}.
This object contain the estimated parameters of a time series process.

}


\section{Slots}{
  \describe{
  
    \item{\code{call}:}{Object of class \code{"call"}: 
        the call of the \code{gsFit} function.
        } 
        
    \item{\code{formula}:}{Object of class \code{"formula"}: 
        a formula object specifying mean and variance      equation.
        }
        
    \item{\code{method}:}{Object of class \code{"character"}: 
        a string denoting the optimization method.
        }
    
    \item{\code{convergence}:}{Object of class \code{"numeric"}:
        an integer code. 0 indicates successful convergence of 
        the of the estimation routine. This variable is the same 
        returned by the one of the optimization methods, which could 
        be one of \code{solnp} or {nlminb} routines. Notice that sometimes
        the optimization algorithm will return a "true" convergence, even 
        when the optimized negative log-likelihood equals to 1e99. In this 
        case, the variable convergence will be equal to 1 to indicate that 
        convergence was not achieved.
        }
        
    \item{\code{messages}:}{Object of class \code{"list"}:
        a character string giving any additional information collected
        during the estimation of the model.
        }
    
   \item{\code{data}:}{Object of class \code{"numeric"}: 
        a numeric vector containing the data of
        the time series to be estimated.
        }
    
        \item{\code{fit}:}{Object of class \code{"list"}: 
        a list with the results from the parameter estimation.
        }
    
        \item{\code{residuals}:}{Object of class \code{"numeric"}: 
        a numeric vector with the residual values.
        }
    
    \item{\code{h.t}:}{Object of class \code{"numeric"}: 
        a numeric vector with the conditional variances.
        }
        
    \item{\code{sigma.t}:}{Object of class \code{"numeric"}: 
        a numeric vector with the conditional standard deviations.
        }
        
    \item{\code{title}:}{Object of class \code{"character"}: 
        a string with the title.
        }
        
    \item{\code{description}:}{Object of class \code{"character"}: 
        a string with a description.
        }

  }
}
\section{Methods}{

  \describe{
  
      \item{show}{\code{signature(object = "fGEVSTABLEGARCH")}: 
        prints an object of class 'fGEVSTABLEGARCH'.
        }  
	 }
}


\author{

Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano 
and Silvia Regina Costa Lopes.

}


\seealso{

\code{gsFit}

}


\keyword{classes}












