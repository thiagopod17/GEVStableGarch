pkgname <- "GEVStableGarch"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('GEVStableGarch')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("dist-gat")
### * dist-gat

flush(stderr()); flush(stdout())

### Name: gat
### Title: Generalized Asymmetric t Distribution
### Aliases: gat dgat pgat qgat rgat
### Keywords: distribution

### ** Examples


# Simulate Random Values and compare with
# the empirical density and probability functions
# Note: This example was addapted from "sstd {fGarch} R Documentation"

# Configure plot and generate random values
par(mfrow = c(2, 2))
set.seed(1000)
r = rgat(n = 1000)
plot(r, type = "l", main = "GAt Random Values", col = "steelblue")

# Plot empirical density and compare with true density:
hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue")
box()
x = seq(min(r), max(r), length = 201)
lines(x, dgat(x), lwd = 2)

# Plot density function and compare with true df:
plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue",
     ylab = "Probability")
lines(x, pgat(x), lwd = 2)

# Compute quantiles:
# Here we compute the quantiles corresponding to the probability points from 
# -10 to 10 and expect to obtain the same input sequence
round(qgat(pgat(q = seq(-10, 10, by = 0.5))), digits = 6)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("dist-skstd")
### * dist-skstd

flush(stderr()); flush(stdout())

### Name: skstd
### Title: Skew Student's t Distribtuion from Fernandez and Steel (1997)
### Aliases: skstd dskstd pskstd qskstd rskstd
### Keywords: distribution

### ** Examples


# Simulate Random Values and compare with
# the empirical density and probability functions
# Note: This example was addapted from "sstd {fGarch} R Documentation"

# Configure plot and generate random values
par(mfrow = c(2, 2))
set.seed(1000)
r = rskstd(n = 1000)
plot(r, type = "l", main = "Skew Student's t Random Values", col = "steelblue")

# Plot empirical density and compare with true density:
hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue")
box()
x = seq(min(r), max(r), length = 201)
lines(x, dskstd(x), lwd = 2)

# Plot density function and compare with true df:
plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue",
     ylab = "Probability")
lines(x, pskstd(x), lwd = 2)

# Compute quantiles:
# Here we compute the quantiles corresponding to the probability points from 
# -10 to 10 and expect to obtain the same input sequence
round(qskstd(pskstd(q = seq(-10, 10, by = 0.5))), digits = 6)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
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
gev.model = gsFit(data = x , formula = ~arma(1,1)+garch(1,1), cond.dist = "gev")



cleanEx()
nameEx("gsMomentAparch")
### * gsMomentAparch

flush(stderr()); flush(stdout())

### Name: gsMomentAparch
### Title: Computation of moments for several conditional distribution
### Aliases: gsMomentAparch
### Keywords: aparch moments

### ** Examples


# Computation of the Moment E( |Z| - gamma Z) ^ delta for several distributions

gsMomentAparch(cond.dist = "stableS1", shape = 1.1, skew = 0, delta = 1.01, gm = 0.99999)

gsMomentAparch(cond.dist = "gev", shape = -4, skew = 0, delta = 1.4, gm = 0)

gsMomentAparch(cond.dist = "gat", shape = c(1.9,2.3), skew = 0.5, delta = 0.4, gm = 0)

gsMomentAparch(cond.dist = "norm", shape = c(1.9,2.3), skew =1, delta = 11.4, gm = -0.999)

gsMomentAparch(cond.dist = "std", shape = 2.001, skew = -0.5, delta = 2, gm = -0.99)

gsMomentAparch(cond.dist = "sstd", shape = 2.001, skew = 0.11, delta = 2, gm = -0.99)

gsMomentAparch(cond.dist = "skstd", shape = 5.001, skew = 0.11, delta = 3, gm = -0.5)

gsMomentAparch(cond.dist = "ged", shape = 6, skew = 0.11, delta = 5.11, gm = -0.5)




cleanEx()
nameEx("gsSelect")
### * gsSelect

flush(stderr()); flush(stdout())

### Name: gsSelect
### Title: Selects the best model according to goodness-of-fit criteria
### Aliases: gsSelect

### ** Examples


# Best ARMA-GARCH model within the range ARMA(0,0)-GARCH(1,0) to ARMA(0,0)-GARCH(1,1)
# using the Corrected Akaike Information Criteria (AICc)
library(fGarch)
data(dem2gbp)
x = dem2gbp[,1]

model = gsSelect (data = x, order.max = c(0,0,1,1), is.aparch = FALSE, 
          algorithm = "sqp", cond.dist = "gev", selection.criteria = "AIC")




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
### Title: Specification of ARMA-GARCH/APARCH models with GEV or stable
###   distributions
### Aliases: gsSpec
### Keywords: models

### ** Examples


# stable-GARCH from Curto et al. (2009) for the DJIA dataset
spec.stable = gsSpec(model = list(mu = 0.0596, omega = 0.0061, 
alpha = 0.0497, beta = 0.9325, skew = -0.9516, shape = 1.9252), 
cond.dist = "stableS1")
sim.stable = gsSim(spec = spec.stable, n = 1000)
 
# GEV-GARCH model from Zhao et al. (2011)
spec.gev = gsSpec(model = list(mu = 0.21, a = 0.32, omega = 0.01,
alpha = 0.45, beta = 0.08, shape = 0.08), cond.dist = "gev")
sim.gev = gsSim(spec = spec.gev, n = 1000)




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
