pkgname <- "GEVStableGarch"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('GEVStableGarch')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("fGEVSTABLEGARCH-class")
### * fGEVSTABLEGARCH-class

flush(stderr()); flush(stdout())

### Name: fGEVSTABLEGARCH-class
### Title: Class '"fGEVSTABLEGARCH"'
### Aliases: fGEVSTABLEGARCH-class show,fGEVSTABLEGARCH-method
### Keywords: classes

### ** Examples

showClass("fGEVSTABLEGARCH")



cleanEx()
nameEx("fGEVSTABLEGARCHSPEC-class")
### * fGEVSTABLEGARCHSPEC-class

flush(stderr()); flush(stdout())

### Name: fGEVSTABLEGARCHSPEC-class
### Title: Class '"fGEVSTABLEGARCHSPEC"'
### Aliases: fGEVSTABLEGARCHSPEC-class show,fGEVSTABLEGARCHSPEC-method
### Keywords: classes

### ** Examples

showClass("fGEVSTABLEGARCHSPEC")



cleanEx()
nameEx("gsFit")
### * gsFit

flush(stderr()); flush(stdout())

### Name: gsFit
### Title: Estimation of ARMA-GARCH/APARCH models
### Aliases: gsFit

### ** Examples

# This examples uses the dataset of the package fGarch to estimate
# an ARMA(1,1)-GARCH(1,1) with GEV conditional distribution.
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]
gev.model = gsFit(data = x , formula = ~garch(1,1), cond.dist = "norm")



cleanEx()
nameEx("gsMomentAparch")
### * gsMomentAparch

flush(stderr()); flush(stdout())

### Name: gsMomentAparch
### Title: Evaluate the moments expression E( |Z| - gamma * Z) ^ delta
### Aliases: gsMomentAparch
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cond.dist = c("stableS1", "gev", "GAt", "norm", "std", 
    "sstd", "skstd", "ged"), shape = 1.5, skew = 0, delta = 1, 
    gm = 0) 
{
    cond.dist = match.arg(cond.dist)
    if (cond.dist == "stableS1") 
        kappa = .stableS1MomentAparch(shape = shape, skew = skew, 
            delta = delta, gm = gm)
    if (cond.dist == "gev") 
        kappa = .gevMomentAparch(shape = shape, delta = delta, 
            gm = gm)
    if (cond.dist == "GAt") 
        kappa = .GAtMomentAparch(shape = shape, delta = delta, 
            skew = skew, gm = gm)
    if (cond.dist == "norm") 
        kappa = .normMomentAparch(delta = delta, gm = gm)
    if (cond.dist == "std") 
        kappa = .stdMomentAparch(shape = shape, delta = delta, 
            gm = gm)
    if (cond.dist == "sstd") 
        kappa = .sstdMomentAparch(shape = shape, skew = skew, 
            delta = delta, gm = gm)
    if (cond.dist == "skstd") 
        kappa = .skstdMomentAparch(shape = shape, skew = skew, 
            delta = delta, gm = gm)
    if (cond.dist == "ged") 
        kappa = .gedMomentAparch(shape = shape, delta = delta, 
            gm = gm)
    kappa
  }



cleanEx()
nameEx("gsSelect")
### * gsSelect

flush(stderr()); flush(stdout())

### Name: gsSelect
### Title: Selects the best model according to goodness-of-fit
### Aliases: gsSelect

### ** Examples

# AIC fit using models from ARMA(0,0)-GARCH(1,0) to ARMA(1,1)-GARCH(1,1)
# with GEV conditional distribution
#library(fGarch)
#data(dem2gbp)
#x = dem2gbp[, 1]
# GSgarch.FitAIC(data = x,1,0,1,0,cond.dist = "gev")



cleanEx()
nameEx("gsSim")
### * gsSim

flush(stderr()); flush(stdout())

### Name: gsSim
### Title: Simulation of ARMA-GARCH/APARCH process
### Aliases: gsSim

### ** Examples

# Simulation of a ARMA-APARCH process with stable conditional distribution
#x <- GSgarch.Sim(N = 2500, mu = 0.1,a = c(0.2,0.3),b = c(0.2,0.5),
#omega = 0.1, alpha = c(0.1,0.2),beta = c(0.1,0.1),gm=c(0.3,-0.3),
#delta = 1,skew = 0.3,shape = 1.9, cond.dis = "stable")



cleanEx()
nameEx("gsSpec")
### * gsSpec

flush(stderr()); flush(stdout())

### Name: gsSpec
### Title: Specification of ARMA-GARCH/APARCH models
### Aliases: gsSpec
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
