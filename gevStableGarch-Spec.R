
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:               SPECIFICATION:
#  GSgarchSpec               Creates a 'fGEVSTABLEGARCHSPEC' object from scratch
################################################################################


# GSgarchSpec function
GSgarchSpec <-
    function (model = list(), presample = NULL,
    cond.dist = c("stable", "gev", "t3", "norm", "std", "sstd", "skstd", "ged"),
    rseed = NULL)
{	
    # Description:
    #   Creates a "GSgarchSpec" object.

    # Arguments:
    #   model - a list with the model parameters as entries
    #     omega - the variance value for GARCH/APARCH
    #       specification,
    #     alpha - a vector of autoregressive coefficients
    #       of length p for the GARCH/APARCH specification,
    #     gm - a vector of leverage coefficients of
    #       length p for the APARCH specification,
    #     beta - a vector of moving average coefficients of
    #       length q for the GARCH/APARCH specification,
    #     mu - the mean value for ARMA specification,
    #     ar - a vector of autoregressive coefficients of
    #       length m for the ARMA specification,
    #     ma - a vector of moving average coefficients of
    #       length n for the ARMA specification,
    #     delta - the exponent value used in the variance equation.
    #     skew - a numeric value listing the distributional
    #        skewness parameter.
    #     shape - a numeric value listing the distributional
    #        shape parameter.
    #   presample - a numeric "matrix" with 3 columns and
    #       at least max(m,n,p,q) rows. 
    #       The first culumn are the innovations, the second
    #       the conditional variances, and the last the time series.
    #       When presample is missing, we construct our presample matrix
    #       as [z,h,y] where z = rnorm(0,1), h = "uev" recursion 
    #       initialization described in Wuertz et al. (2006) and
    #       y = mu. Note that the conditional variance column
    #       can contain only strictly positive numbers.
    #       If the model is a pure AR or MA the presample 
    #       matrix will become a simple array.
    #   cond.dist - a character string naming the distribution
    #       function.
    #   rseed - optional random seed. The default seed is '0'.
     
    # Return: An object of class fGEVSTABLEGARCHSPEC   
        # Slots:
        #   call - the function call.
        #   formula - a formula object describing the model, e.g.
        #       ARMA(m,n) + GARCH(p,q). ARMA can be missing or
        #       specified as AR(m) or MA(n) in the case of pure
        #       autoregressive or moving average models. GARCH may
        #       alternatively be specified as ARCH(p) or APARCH(p,q).
        #   model - as declared in the input.

    # FUNCTION
      
    # Skewness Parameter Settings:
    skew = list(
        "dstable" = 0,
        "dgev" = NULL,
        "dt3" = 1,
        "dnorm" = NULL,
        "dstd" = NULL,
        "dsstd" = 0.9,
        "dskstd" = 1,
        "dged" = NULL)

    # Shape Parameter Settings:
    shape = list(
      "dstable" = 1.7,
      "dgev" = 0.3,
      "dt3" = c(3,1),
      "dnorm" = NULL,
      "dstd" = 4,
      "dsstd" = 4,
      "dskstd" = 4,
      "dged" = 4)

    # Conditional distribution
    cond.dist = match.arg(cond.dist)

    # Default Model (AR(1)):
    initialDelta = NULL
    if(!is.null(model$alpha) && is.null(model$delta)) # Garch model
    {
        if(cond.dist == "dstable")
            initialDelta = 1
        else 
            initialDelta = 2
    }
    	
    control = list(
        omega = 1,
        alpha = NULL,
        gm = NULL,
        beta = NULL,
        mu = NULL,
        ar = NULL,
        ma = NULL,
        delta = initialDelta,
        skew = skew[[cond.dist]],
        shape = shape[[cond.dist]]
        )

    # Update Control:
    control[names(model)] <- model
    model <- control

    # Model Orders:
    order.ar = length(model$ar)
    order.ma = length(model$ma)
    order.alpha = length(model$alpha)
    order.beta = length(model$beta)
    order.max = max(order.ar, order.ma, order.alpha, order.beta)
    
    
    # Error treatment of input parameters:   
    if((order.alpha == 0 && order.beta != 0))
    	stop("In a Garch(p,q)/Aparch(p,q) model we must have p > 0")  
       
    if(!is.null(model$delta)){
    	if(length(model$delta) > 1)
       stop("The parameter 'delta' must be a single number")
    	 if(!(model$delta > 0))
      	stop("The parameter 'delta' must be > 0.")
    }
    
    if( (length(model$gm) != 0) && (length(model$alpha) != 0)) # means aparch model
    {
      if(length(model$alpha) != length(model$gm))
        stop("'alpha' and 'gm' must have the same size for APARCH models")    
    }
   
    if(!is.null(model$gm)){
    	if(sum(!(abs(model$gm)<1)) > 0) # all gm in (-1,1)
        stop("The parameter 'gm' must be in the range -1 < gm < 1")     	
    }
             
    if(sum(model$alpha < 0) > 0) # all alpha in [0,+infty)
        stop("The parameter 'alpha' must be greater than or equal 0")          
        
        
    if(sum(model$beta < 0) > 0) # all beta in [0,+infty)
        stop("The parameter 'beta' must be greater than or equal 0")
        
    if(!(model$omega > 0)) # omega > 0
        stop("The parameter 'omega' must be > 0")        
            
    
    if(is.null(model$ar) && is.null(model$ma) && is.null(model$alpha))
    	stop("The model parameters were not specified correctly.")
    	
    if(!is.null(model$delta) && is.null(model$alpha)) 	
       stop("The parameter delta should be only specified for a GARCH or APARCH model.")
       
        # if we have a presample check if it has the correct range.       
    if(is.matrix(presample) && !is.null(presample) && !is.null(model$alpha) ){
    	if((dim(presample)[1] != order.max) || (dim(presample)[2] != 3) )
    	   stop("The presample object should be a matrix with three columns representing the innovations, conditional variance and the time series, respectively. The program expects a matrix of dimensions l x 3, where l = max(m,n,p,q)")
    } 
    if(is.matrix(presample) && !is.null(presample) && is.null(model$alpha) ){
    	if((dim(presample)[1] != order.max) || (dim(presample)[2] != 2) )
    	   stop("The presample object should be an array with the time series. The program expects an array with dimensions l x 2, where l = max(m,n,p,q)")
    }       
       
       	
    # Compose Mean Formula Object:
    formula.mean = ""
    if (order.ar == 0 && order.ma == 0) {
        formula.mean = ""
    }
    else {
        formula.mean = paste ("arma(", as.character(order.ar), ", ",
            as.character(order.ma), ")", sep = "")
    }
    
    
#     if (order.ar > 0 && order.ma == 0) {
#         formula.mean = paste ("ar(", as.character(order.ar), ")", sep = "")
#     }
#     if (order.ar == 0 && order.ma > 0) {
#         formula.mean = paste ("ma(", as.character(order.ma), ")",  sep = "")
#     }
#     if (order.ar > 0 && order.ma > 0) {
#         formula.mean = paste ("arma(", as.character(order.ar), ", ",
#             as.character(order.ma), ")", sep = "")
#     }

    # Compose Variance Formula Object:
    formula.var = ""
    #if (order.alpha > 0 && order.beta == 0) formula.var = "arch"
    if (order.alpha > 0 && order.beta == 0) formula.var = "garch"
    if (order.alpha > 0 && order.beta > 0) formula.var = "garch"   
    if(!is.null(model$alpha)){
    	if (!is.null(model$gm)){
    		if(sum(model$gm == 0) != length(model$gm))
    		   formula.var = "aparch"
    	}
    	if (model$delta != 2 && cond.dist != "stable") 
          formula.var = "aparch" # gm = 0 and delta != 0 we get powergarch model   	
    	if (model$delta != 1 && cond.dist == "stable") 
    	    formula.var = "aparch" 
    }
   
    if (order.alpha == 0 && order.beta == 0) {
        formula.var = formula.var
    }
#     if (order.alpha > 0 && order.beta == 0) {
#         formula.var = paste(formula.var, "(", as.character(order.alpha), ")",
#             sep = "")
#     } 
#     if (order.alpha > 0 && order.beta > 0) {
#         formula.var = paste(formula.var, "(", as.character(order.alpha),
#             ", ", as.character(order.beta), ")", sep = "")
#     }
      if (order.alpha > 0) {
        formula.var = paste(formula.var, "(", as.character(order.alpha),
                            ", ", as.character(order.beta), ")", sep = "")
      }

    # Compose Mean-Variance Formula Object:
    if (formula.mean == "") {
        formula = as.formula(paste("~", formula.var))
    } 
    if (formula.var == "") {
        formula = as.formula(paste("~", formula.mean))
    }
    if ((formula.mean != "") && (formula.var != "")){
        formula = as.formula(paste("~", formula.mean, "+", formula.var))
    }

    # Add NULL default entries:
    if (is.null(model$mu)) model$mu = 0
    if (is.null(model$ar)) model$ar = 0
    if (is.null(model$ma)) model$ma = 0
    if (is.null(model$gm)) model$gm = rep(0, times = order.alpha)

    # Seed:
    if (is.null(rseed)) {
        rseed = 0
    } else {
        set.seed(rseed)
    }

    # Define Missing Presample:
    if(formula.var != ""){
	 persistency = 1-sum(model$alpha)-sum(model$beta)
	 if(persistency*(1-persistency) < 0) # avoid to construct a presample with negative conditional variance.
	    persistency = 0.1 
	 if(!is.null(presample) && is.matrix(presample)){
	     z = presample[, 1]
	     h = presample[, 2]
	     y = presample[, 3]
	     if(sum(!(h > 0)))
	        stop("Conditional Variance column can have only strictly positive numbers")
	 }else{
	     z = rnorm(n = order.max)
	     h = rep(model$omega*(1 + persistency*(1-persistency)), times = order.max)
	     y = rep(model$mu, times = order.max)
	 }
	 presample = cbind(z, h, y)
    }else{	
       if(!is.null(presample))
       {
	     z = presample[, 1]
	     y = presample[, 2]          
       }
       else
       {
       	 z = rnorm(n = order.max)
         y = rep(0,order.max)
       }
       presample = cbind(z,y)
    }
    
    # Result: 
    new("fGEVSTABLEGARCHSPEC",
        call = match.call(),
        formula = formula,
        model = list(omega = model$omega, alpha = model$alpha,
            gm = model$gm, beta = model$beta, mu = model$mu,
            ar = model$ar, ma = model$ma, delta = model$delta,
            skew = model$skew, shape = model$shape),
        presample = as.matrix(presample),
        distribution = as.character(cond.dist),
        rseed = as.numeric(rseed)
    )
}


################################################################################

