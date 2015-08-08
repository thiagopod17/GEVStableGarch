
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
# TEST THE MAIN FUNCTIONS OF PACKAGE GEVStableGarch: 
# 
# Note: There is also the possibility to enforce the tests when 
# checking the package with R CMD.
################################################################################

# Required packages
library(testthat)
library(fGarch)
library(GEVStableGarch)
library(Rsolnp)
library(stable)
library(stabledist)
library(skewt)

# Create file to save the result of the test
log.file.directory = '/Users/thiago/dropbox/Pappers/GEVStableGarch/GEVStableGarch/GEVStableGarch-unitTest-results'
log.file.name = format(Sys.time(), "DATE-%m-%d-%y---TIME-%Hh%Mm%Ssec")
log.file.adress = paste(log.file.directory, "/", log.file.name, sep = '')

write.log.function = function(x) {
  write(x, file = log.file.adress,
        ncolumns = if(is.character(x)) 1 else 5,
        append = TRUE, sep = "")
}
write.log.function("===================================================")
write.log.function("Testing functions from package GEVStableGarch.")
write.log.function(paste("\nDate and time: ",log.file.name))
write.log.function("===================================================")





# Comparing the estimated parameters of a 
# Garch(1,1) with normal distribution and the fGarch package
# Compare the results obtained from 3 different algorithms 
# againts the reported results from garchFit, G@RCH and Splus FinMetrics softwares.

data(dem2gbp)
x = dem2gbp[,1]

garch11.garchFit = c(-0.006190314, 0.010761384, 0.153134059, 0.805973748)
garch11.garch.software = c(-0.006184, 0.010760, 0.153407, 0.805874)
garch11.finmetrics = c(-0.006053, 0.010896, 0.15421, 0.80445)

garch11.gevstablegarch.sqp = gsFit(data = x, formula = ~garch(1,1),
                          cond.dist = "norm", algorithm = "sqp")@fit$par

garch11.gevstablegarch.nlminb = gsFit(data = x, formula = ~garch(1,1),
                                         cond.dist = "norm", algorithm = "nlminb")@fit$par

garch11.gevstablegarch.nlminb.nm = gsFit(data = x, formula = ~garch(1,1),
                                   cond.dist = "norm", algorithm = "nlminb+nm")@fit$par

garch11.gevstablegarch.sqp.restriction = gsFit(data = x, formula = ~garch(1,1),
                                            cond.dist = "norm", algorithm = "sqp.restriction")@fit$par

write.log.function("\n\n-----\nTesting the estimated parameters from the 
garch(1,1)-normal against garchFit, G@RCH and SPLUS Fin Metrics...")
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.garchFit)/garch11.gevstablegarch.sqp), 1e-4)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.garch.software)/garch11.gevstablegarch.sqp), 1e-2)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.garchFit)/garch11.gevstablegarch.sqp), 1e-4)
write.log.function("OK")


write.log.function("\n\n-----\nTesting he estimated parameters from the garch(1,1)-normal 
algorithms: sqp, nlminb, nlminb+nm and sqp.restriction...")
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.gevstablegarch.nlminb)/garch11.gevstablegarch.nlminb), 1e-4)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.gevstablegarch.nlminb)/garch11.gevstablegarch.nlminb.nm), 1e-4)
expect_less_than( max(abs(garch11.gevstablegarch.nlminb - garch11.gevstablegarch.nlminb)/garch11.gevstablegarch.nlminb.nm), 1e-4)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.gevstablegarch.sqp.restriction)/garch11.gevstablegarch.sqp), 1e-4)
write.log.function("OK")







