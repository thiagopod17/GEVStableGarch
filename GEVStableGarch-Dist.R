
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
#  GSgarch.Dist        Computes density values for several 
#                      conditional distributions
################################################################################


GSgarch.Dist <-
    function(z, hh, shape = 4, skew = 0.1, 
    cond.dist = "sstd", GStol = 1e-8) 
{
    # Description:
      
      ### Constants to adjust tolerances ###
      # GStol <- 1e-8 # general tolerance for arma-garch parameters. In the beggining it was set to 1e-5
      # GStol.b <- 1e-7 # upper and lower bounds tolerance. Should be greater than tol
      # GSstable.tol <- 1e-2 # tolerance for stable distribution parameter set
      # GSstable.tol.b <- 2e-2 # boundary tolerance. Should be greater than GSstable.tol
      
      
      
      ###########################
      ###  Time series model  ###
      # According to Wurtz et al. (2006) 
      # xt = ? + a(B)xt + b(B)et,
      # et = zt * sigmat
      # zt ~ Dv(0, 1)
      # sigmat^2 = omega + sum(alphai(e(t-i))^2) + sum(betaj*sigmat-j^2)
      
      #################################################
      ###  Variables notation used inside functions ###
      # N: sample size
      # m,n: ARMA order
      # p,q, pq = max(p,q): GARCH order. They were called u,v, uv in Wurtz et al.(2006)
      # mu, a, b: ARMA parameters
      # omega, alpha, gamma, beta: APARCH vector parameters
      # garchLLH: Log Likelihood for the ARMA(m,n)-GARCH(p,q) model
      # parm: ARMA-GARCH parameters concatenated as c(mu,a,b,omega,alpha,gamma,beta)
      # h: conditional variance. It was referenced as sigmat in model equations
      # x: data set. It is the time series that will be modelled as ARMA-GARCH
      # xi: GEV shape parameter. xi > -0.5 (llh) and < 0.5 (finiteness of variance, see Zhao et al. (2011))
      # AR, MA: if AR (or MA) equals TRUE, then the model has AR (or MA) order equals to zero.
      # param or pm (stabledist): Parametrization of stable distributions as chosen to be 2.  
      # Even for stable distribution llh, there's a problem for finding the estimators for
      # alpha near to 2 and beta > 0. The bad performance on the ARMA-GARCH model was
      # similar to the single llh estimation of i.i.d samples of stable distribution. 
      # In our simulations, the param = 2 was chosen to be the ones that performs
      # better in these situations.       
      
      #################
      #### Remarks ####
      # et for ARMA(m,n), n > 1 will be initiated as 0, i.e, e[-n+1:0] = 0.
      # ht will be initiated as 0.1 (see eq. (22) of Wurtz, 2006).
      
    #   We choose not to avoid calling this function with out of bound parameters
    #   because in these cases we simply return a huge likelihood to guide the optimization algorithm.
    cond.dist.list <- c("norm", "std", "sstd", "gev", "stable")
    if( !any(cond.dist.list == cond.dist) )   
      stop ("Invalid Conditional Distribution. Choose: norm,std,sstd,gev or stable")
    if(sum(is.na(hh)) > 0 || min(hh) == 0)
    {
      #stop ("NA or zero element found in vector hh")
      warning("NA or zero element found in vector hh")
      return(1e99)
    }
    
    # Processing
    if(cond.dist == "norm")
      return(-sum(log(dnorm(x = z/hh)/hh)))
    if(cond.dist == "std")
    {
      if(!(shape > 2))
        stop("Invalid shape in std distribution. shape > 2")
      nu = shape
      return(-sum(log(dstd(x = z/hh, nu = nu)/hh)))
    }
    if(cond.dist == "sstd")
    {
      if(!(shape > 2) || !(skew > 0))
        stop("Invalid shape or skew in skewt parameters. shape > 2 and skew > 0")
      gm = skew
      xi = gm
      nu = shape
      return(-sum(log(dsstd(x = z/hh, nu = nu, xi = xi)/hh)))
    }
    if(cond.dist == "gev")
    {
      if(abs(shape) < 1e-6)
        stop("shape parameter from GEV is to small. abs(shape) < 1e-6")
      sig <- hh
      xi <- shape
      arg <- z/sig
      y <- 1 + xi * arg
      if(sum(is.na(y)) || sum(is.nan(y)) ||
           sum(is.infinite(y))) 
        return(1e99)
      gev.cond.gD <- any(y < GStol ) || abs(xi) < 1e-6
      if(gev.cond.gD)
      {
        return(1e99)
      }
      llh <- sum(log(sig)) + sum(y^(-1/xi)) + sum(log(y))*(1/xi + 1)
      return(llh)
    }
    if(cond.dist == "stable")
    {
      if( !(shape > 1) || !(shape < 2) || !(abs(skew) < 1))
        stop("Invalid shape or skew in stable parameters. 1 < shape < 2 and |skew| < 1")
      sig <- hh
      alpha.stable <- shape
      beta.stable <- skew
      arg <- z/sig
      y <- arg
      if(sum(is.na(y)) || sum(is.nan(y)) ||
           sum(is.infinite(y)))
        return(1e99)
      #dens.stable <- stable::dstable.quick(y,alpha = alpha.stable, 
      #beta = beta.stable, gamma = 1, delta = 0, param = 2)
      dens.stable <- GSgarch.dstable(y,alpha = alpha.stable,
                                     beta = beta.stable, gamma = 1,delta = 0, param = 2)
      llh <- sum(log(sig[sig>0])) - 
        sum(log(dens.stable[dens.stable>0]))
      return(llh)
    }
}




################################################################################
