
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
# TEST CASES FOF FUNCTION:               SPECIFICATION:
#  GSgarch.Fit               Useful for finding erros. Please add the expected 
#  						               output of the test case whenever possible. 
#							               This will help to trace erros during debugging.
################################################################################


# We will make the interface of this function more similar to the garchFit function
# from package fGarch to make its use easier. 
# The input parameters for this version (march/2015) are:

  # OK. algorithm: a string parameter that determines the algorithm used for maximum likelihood estimation.
  # OK. cond.dist: name of the conditional distribution, one of gev, stable, norm, std, sstd
  # OK. control: control parameters, the same as used for the functions from nlminb, and 'bfgs' and 'Nelder-Mead' from optim.
  # OK. data: The dataset to be estimated.
  # formula: formula object describing the mean and variance equation of the ARMA-GARCH/APARCH model.
  # OK. intercept: this flag determines if the parameter for the mean will be estimated or not
  # OK. print.Result (Padrao eh TRUE): A boolean variable specifying whether or not the user wants to print the results after the function calling.
  # OK. get.res: (NAO VAMOS TER MAIS ESSA VARIAVEL)
  # OK. GSstable.tol e GStol: (CONFIGURAR NO INICIO DA FUNCAO, NAO MAIS NECESSARIA AQUI)
  # APARCH: Vamos tirar pois usaremos da formula.


?garchFit

library(GEVStableGarch)
?GSgarch.Fit

# General Testing of invalid input parameters with cond.dist, algorithm and data (with NA,
# inf, NULL or character)
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]
gF.new = GSgarch.Fit(data = x , m = 1,n = 1,p = 1,q = 1,
                     cond.dist = "gev", intercept = TRUE, APARCH = TRUE, 
                     algorithm = "sqp.restriction")

# Testing different type of datasets as input
x <- c("asdf",rnorm(1000))
x <- c(rnorm(100),NA,rnorm(200))
x <- c(rnorm(100),-Inf)
x <- c(NULL)
x <- rnorm(100)

GSgarch.Fit(data = x , m = 1,n = 1,p = 1,q = 1,
            cond.dist = "gev", intercept = TRUE, APARCH = FALSE, 
            algorithm = "sqp.restriction")

# ARMA(1,1)-GARCH(1,1) with and without intercept 
  data = eval(parse(text = data))
x = dem2gbp[, 1]
gF.new <- GSgarch.Fit(data = x[1:200] , m = 1,n = 1,p = 1,q = 1,
            cond.dist = "gev", intercept = FALSE, APARCH = FALSE, 
            algorithm = "sqp.restriction")

# GARCH(1,1) with and without intercept 
x = dem2gbp[, 1]
gF.new <- GSgarch.Fit(data = x , m = 0,n = 0,p = 1,q = 1,
                      cond.dist = "norm", intercept = FALSE, APARCH = FALSE, 
                      algorithm = "sqp")

# Understanding the 'forumla' objects
formula = arma(2,1)~garch(4,1)
# length 3
length(formula)
formula = ~arma(1,1)
formula = ~ar(1)
formula = ~aparch(1,1)
formula = ~garch(1,1)
formula = ~arch(1)
length(formula)
# length 2

################
# Debugging function '.garchArgsParser' from package fGarch
data = dem2gbp[, 1]
formula = arma(2,1)~garch(4,1)

# Commands before calling this function
# Parse formula and data for garchFit ...
#   Note in the new version we are working with timeSeries ...

TESTE1 <- function(formula, data)
{
    Name = capture.output(substitute(data))
    if(is.character(data)) {
      eval(parse(text = paste("data(", data, ")")))
    }
    # data <- if (inherits(data, "timeSeries") data else as.timeSeries(data)
    data <- as.data.frame(data)
    
    # Column Names:
    if (isUnivariate(data)) {
      colnames(data) <- "data"
    } else {
      # Check unique column Names:
      uniqueNames = unique(sort(colnames(data)))
      if (is.null(colnames(data))) {
        stop("Column names of data are missing.")
      }
      if (length(colnames(data)) != length(uniqueNames)) {
        stop("Column names of data are not unique.")
      }
    }
    
    # Handle if we have no left-hand-side for the formula ...
    #   Note in this case the length of the formula is 2 (else 3):
    if (length(formula) == 3 && isUnivariate(data) ) formula[2] <- NULL
    if (length(formula) == 2) {
      if (isUnivariate(data)) {
        # Missing lhs -- we substitute the data file name as lhs ...
        formula = as.formula(paste("data", paste(formula, collapse = " ")))
      } else {
        stop("Multivariate data inputs require lhs for the formula.")
      }
    }
    garchArgsParser2(formula = formula, data = data, trace = TRUE)
}
garchFit(formula = ~arma(1,1)+garch(1,1))
garchArgsParser2 <- 
?garchFit
function (formula, data, trace = FALSE) 
{
  allVars = unique(sort(all.vars(formula)))
  allVarsTest = mean(allVars %in% colnames(data))
  if (allVarsTest != 1) {
    print(allVars)
    print(colnames(data))
    stop("Formula and data units do not match.")
  }
  formula.lhs = as.character(formula)[2]
  mf = match.call(expand.dots = FALSE)
  if (trace) {
    cat("\nMatched Function Call:\n ")
    print(mf)
  }
  m = match(c("formula", "data"), names(mf), 0)
  mf = mf[c(1, m)]
  mf[[1]] = as.name(".garchModelSeries")
  mf$fake = FALSE
  mf$lhs = TRUE
  if (trace) {
    cat("\nModelSeries Call:\n ")
    print(mf)
  }
  x = eval(mf, parent.frame())
  if (trace) 
    print(x)
  x = as.vector(x[, 1])
  names(x) = rownames(data)
  if (trace) 
    print(x)
  allLabels = attr(terms(formula), "term.labels")
  if (trace) {
    cat("\nAll Term Labels:\n ")
    print(allLabels)
  }
  if (length(allLabels) == 2) {
    formula.mean = as.formula(paste("~", allLabels[1]))
    formula.var = as.formula(paste("~", allLabels[2]))
  }
  else if (length(allLabels) == 1) {
    formula.mean = as.formula("~ arma(0, 0)")
    formula.var = as.formula(paste("~", allLabels[1]))
  }
  if (trace) {
    cat("\nMean Formula:\n ")
    print(formula.mean)
    cat("\nVariance Formula:\n ")
    print(formula.var)
  }
  ans <- list(formula.mean = formula.mean, formula.var = formula.var, 
              formula.lhs = formula.lhs, series = x)
  ans
}

TESTE1(formula = formula, data = data)


##########
allVars = unique(sort(all.vars(formula)))
formula.lhs = as.character(formula)[2]
allLabels = attr(terms(formula), "term.labels")



if (length(allLabels) == 2) {
  formula.mean = as.formula(paste("~", allLabels[1]))
  formula.var = as.formula(paste("~", allLabels[2]))
}
if (length(allLabels) == 1) {
  formula.mean = as.formula("~ arma(0, 0)")
  formula.var = as.formula(paste("~", allLabels[1]))
}
