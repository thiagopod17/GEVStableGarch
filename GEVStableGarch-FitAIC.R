
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
#  FUNCTION:               DESCRIPTION:
#
#  GSgarch.FitAIC          Find best fitted model according to AIC criterion
################################################################################

GSgarch.FitAIC <- 
    function(data,mMAX=1,nMAX=1,pMAX=1,qMAX=1, cond.dist = "norm", 
    algorithm = "sqp",APARCH = FALSE, intercept = TRUE,control = NULL)
{
  # Description:
  #   error treatment on input parameters
  cond.dist.list <- c("norm", "std", "sstd", "gev", "stable")
  if( !any(cond.dist.list == cond.dist) )   
    stop ("Invalid Conditional Distribution. Choose: norm,std,sstd,gev or stable")
  if( !is.numeric(data) || !is.vector(data))
    stop("data set must be a numerical one dimensional vector")
  
  # begin of function
  T <- length(data)
  aic.min <- 1e10
  aic <- 1e10
  fit.min <- list()
  fit <- list()
  aic.list <- GSgarch.GetOrder(mMAX,nMAX,pMAX,qMAX)
  aic.list.size <- length(aic.list[,1])
  for( i in 1:aic.list.size)
  {
    fit = GSgarch.Fit(data, m = aic.list[i,1],n = aic.list[i,2],
                      p = aic.list[i,3],q = aic.list[i,4], cond.dist = cond.dist, 
                      APARCH = APARCH, algorithm = algorithm, 
                      intercept = intercept, control = control)
    nParam <- length(fit$par)
    aic  = 2*fit$llh + 2*nParam 
    aicc = 2*fit$llh + 2*nParam*T/(T - nParam - 1)
    bic =  2*fit$llh + nParam*log(T)
    cat(aic.list[i,],"-llh:",fit$llh,"AIC:",aic," AICC:",aicc," BIC:",bic,"\n")
    if (aic < aic.min)
    {
      fit.min <- fit
      fit.min$order <- aic.list[i,]
      aic.min <- aic
    }
  }
  return(fit.min)
}



################################################################################

