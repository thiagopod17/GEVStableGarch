
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
#  .armaGarchDist          Computes density values for several 
#                          conditional distributions
################################################################################


.armaGarchDist <- 
    function(z, hh, shape = 1.5, skew = 0, 
    cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "GAt", "norm", "std", "sstd", "skstd", "ged"), 
    TOLG = 1e-8) 
{
    # Description:
    #   Calculates the likelihood function for a vector of points (z)
    #   according to the specified distribution.
    #   REMARKS: We choose not to avoid calling this function with out of bound parameters
    #   because in these cases we simply return a huge likelihood to guide 
    #   the optimization algorithm.
    
    # Arguments:
    #   z - vector of points to calculate the llh.
    #   h - vector of conditional variances ????? BETTER DESCRIPTION NEEDED.
    #   cond.dist - name of the conditional distribution
    #   TOLG - general tolerance for arma-garch parameters. 
    #   In the beggining it was set to 1e-5
      
    # FUNCTION:  
      
    # Error treatment of input parameters
    cond.dist = match.arg(cond.dist)
    if(length(z) != length(hh))
        stop("Error: Vectors 'z' and 'hh' have different length.")
    if(sum(is.na(hh)) > 0 || min(hh) == 0)
    {
        return(1e99)
    }
    
    # normal conditional distribution
    if(cond.dist == "norm")
        return(-sum(log(dnorm(x = z/hh)/hh)))
    
    # t-student conditional distribution.
    if(cond.dist == "std")
    {
        if(!(shape > 2))
            return(Inf)
        return(-sum(log(dstd(x = z/hh, nu = shape)/hh)))
    }
    
    # skew t-student conditional (standardized version defined in Wurtz)
    if(cond.dist == "sstd")
    {
        if(!(shape > 2) || !(skew > 0))
            return(1e99)       
        return(-sum(log(dsstd(x = z/hh, nu = shape, xi = skew)/hh)))        
 
#         M1 = sqrt((shape-2)/pi)*gamma(shape/2)^(-1)*
#           gamma((shape-1)/2)
#         M2 = 1
#         return(-sum(log(dsstd(x = z/hh, nu = nu, xi = xi, mean = (skew-1/skew)*M1,
#                sd = sqrt((M2-M1^2)*(skew^2+1/skew^2)+2*M1^2-M2))/hh)))
    }
    
    # skew t-student from Fernandez, C. and Steel, M. F. J. (1998)
    if(cond.dist == "skstd")
    {
      if(!(shape > 2) || !(skew > 0))
          return(1e99)     
      return(-sum(log(dskstd(x = z/hh, nu = shape, xi = skew)/hh)))        
    
    }
    
    # GAt distribution
    if(cond.dist == "GAt")
    {
        if(!(shape[1] > 0) || !(shape[2] > 0) || !(skew > 0))
           return(1e99)    
        return(-sum(log(dGAt(x = z/hh, nu = shape[1], d = shape[2], xi = skew)/hh)))        
    }
    
    # GED conditional distribution.
    if(cond.dist == "ged")
    {
      if(!(shape > 0))
          return(1e99)
      return(-sum(log(dged(x = z/hh, nu = shape)/hh)))
    }
    
    # GEV conditional distribution
    if(cond.dist == "gev")
    {
       if( (shape[1] <= -0.5) || (shape[1] >= 0.5)) # to ensure good mle properties and finiteness of the variance
            return(Inf)
       #if(sum(is.na(z/hh)) || sum(is.nan(z/hh)) || sum(is.infinite(z/hh))) 
       #     return(Inf)
        result = -sum(log( (dgev(x = z/hh, xi = shape[1]) + shape[2])/hh))
        #print("==============")
        #print(result)
        #print(shape)
        #print(sum(z/hh))
        return(result)
#         if(abs(shape) < 1e-6)
#           stop("shape parameter from GEV is to small. abs(shape) < 1e-6")
#        sig <- hh
#        xi <- shape
#        arg <- z/sig
#        y <- 1 + xi * arg
#        if(sum(is.na(y)) || sum(is.nan(y)) || sum(is.infinite(y))) 
#	  {
#            print("sum(is.na(y)) || sum(is.nan(y)) || sum(is.infinite(y))")
#            return(1e99)
#        }
#        gev.cond.gD <- any( y < TOLG ) || abs(xi) < 1e-6
#        if(gev.cond.gD)
#        {
#          print("any( y < TOLG ) || abs(xi) < 1e-6")
#          #return(1e99)
#        }
#        llh <- sum(log(sig)) + sum(y^(-1/xi)) + sum(log(y))*(1/xi + 1)
        # print(llh)
 #       return(llh)
    }
    
    # stable conditional distribution
    if( any ( cond.dist == c("stableS0", "stableS1", "stableS2") ) )
    {
        # Return Big Values if we are out of parameter space

        if( !(shape > 1) || !(shape < 2) || !(abs(skew) < 1))
            return(1e99)        

        y <- z/hh
        if(sum(is.na(y)) || sum(is.nan(y)) || sum(is.infinite(y)))
            return(1e99)
        
        # Compute density either using "stable" or "stabledist" package
        #dens.stable = stable::dstable.quick
        
        if(cond.dist == "stableS0")
            result = -sum(log(stable::dstable.quick(x = z/hh, alpha = shape,
                      beta = skew, param = 0)/hh))
        
        if(cond.dist == "stableS1")
          result = -sum(log(stable::dstable.quick(x = z/hh, alpha = shape,
                                        beta = skew, param = 1)/hh))
        
        if(cond.dist == "stableS2")
          result = -sum(log(stable::dstable.quick(x = z/hh, alpha = shape,
                                        beta = skew, param = 2)/hh))
        
        # Result 
        return(result)
        
        
#         dens.stable <- stable::dstable.quick(y,alpha = shape, 
#         beta = skew, gamma = 1, delta = 0, param = 1)
#         #dens.stable <- GSgarch.dstable(y,alpha = shape,
#         #               beta = skew, gamma = 1,delta = 0, param = 2)
#         llh <- sum(log(sig[sig>0])) - 
#             sum(log(dens.stable[dens.stable>0]))
#         return(llh)
    }
}



################################################################################
