
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
#  .getStart           Test on input parameters for different models 
#    					                including pure ARMA models							               
################################################################################


#######
# Test the output for pure arma models
library(fGarch)
data(dem2gbp)
x = dem2gbp[,1]
.getStart(data = x,m = 1,n = 1,p = 0,q = 0, 
                 AR = FALSE, MA = FALSE,
                 cond.dist = "norm", GSstable.tol.b = 2e-2, GStol.b = 1e-7)  


# Test the output for the t3 distribution
library(fGarch)
data(dem2gbp)
x = dem2gbp[,1]
.getStart(data = x,m = 3,n = 3,p = 3,q = 2, 
                 AR = FALSE, MA = FALSE, ARMAonly = TRUE,
                 cond.dist = "t3", GSstable.tol.b = 2e-2, GStol.b = 1e-7)  





