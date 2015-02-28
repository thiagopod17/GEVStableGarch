############



library(fGarch)
data(dem2gbp)
library(GEVStableGarch)
x = dem2gbp[, 1]
res <- garch11Fit(x, start.h = var(x))
res2 <- garchFit(data = x, formula = ~garch(1,1))
res6 <- GSgarch.Fit(data = x , formula = ~garch(1,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "sqp", DEBUG = TRUE)
res8 <- GSgarch.Fit(data = x , formula = ~garch(1,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb", DEBUG = TRUE)


depois <- res2@fit$matcoef[,1]-res8$matcoef[,1]
antes <- res2@fit$matcoef[,1]-res6$matcoef[,1]

res7 <- GSgarch.Fit(data = x, 0,0,1,1, cond.dist = "norm", 
            intercept = TRUE, algorithm = "nlminb")
error2 <- res2@fit$matcoef[,1]-res
error3 <- res2@fit$matcoef[,1]-res3$matcoef[,1]
error4 <- res2@fit$matcoef[,1]-res4$matcoef[,1] # Sem filtro arma, com filtro do Wuertz
error5 # colocando filtro da funcao garch11Fit

depois <- res2@fit$matcoef[,1]-res6$matcoef[,1]
antes <- res2@fit$matcoef[,1]-res7$par
# GSgarch.Fit(formula = ~garch(1,1),data = x)
var(x)/0.2229

res4$matcoef[,1]-res5$matcoef[,1]


# myFilterFunction
m = 1
n = 1
p = 1
q = 1
mn = max(m,n)
pq = max(p,q)
data = x
mu = mean(data)
a = 0
b = 0
omega = 0.1
alpha = 0.1
beta = 0.8
gm = 0
delta = 2
N = length(data)



# # tapper off the mean equation 2
# e.init <- rep(0,mn)
# e.parc <- c(e.init,filter(data, filter = c(1, -a), sides = 1)[(mn+1):N])
# e.res <- c( e.init, filter(e.parc[-(1:mn)], filter = -b,
#                            method = "recursive", init = e.init[1:n]))     
# 
# # find i.i.d sequence z
#   z <- e.res - mu
# 
#   h <- rep(0.1, pq)
#   edeltat = 0
#   for( i in 1:p)
#   {
#     edelta <- alpha[i]*(abs(z)-gm[i]*z)^delta
#     edeltat = edeltat +  edelta[(p-(i-1)):(N-i)]
#   }
#   edeltat = c(h[1:p],edeltat)
#   c <- omega/(1-sum(beta))
#   h <- c( h[1:pq], c + filter(edeltat[-(1:pq)], filter = beta,
#                               method = "recursive", init = h[q:1]-c))
#   hh <- abs(h)^(1/delta)
# 
# e.res
# h
# hh


garch11Fit = function(x, start.h = 0.1)
{
  flag = TRUE
  # Step 1: Initialize Time Series Globally:
  x <<- x
  # Step 2: Initialize Model Parameters and Bounds:
  Mean = mean(x); Var = var(x); S = 1e-6
  params = c(mu = Mean, omega = 0.1, alpha = 0.1, beta = 0.8)
  lowerBounds = c(mu = -10*abs(Mean), omega = S^2, alpha = S, beta = S)
  upperBounds = c(mu = 10*abs(Mean), omega = 100*Var, alpha = 1-S, beta = 1-S)
  
  # Step 3: Set Conditional Distribution Function:
  garchDist = function(z, hh) { dnorm(x = z/hh)/hh }
  
  # Step 4: Compose log-Likelihood Function:
  garchLLH = function(parm) {
    mu = parm[1]; omega = parm[2]; alpha = parm[3]; beta = parm[4]
    z = (x-mu); Mean = mean(z^2)
    
    # Use Filter Representation:
    e = omega + alpha * c(Mean, z[-length(x)]^2)
    #h = filter(e, beta, "r", init = Mean)
    
    
    # Conditional Variance filtering - filter - Wuertz et al. (2006)
    gamma = 0
    delta  = 2
    eps = z
    uv = 1
    h <- rep(start.h, uv)
    edelta = (abs(eps)-gamma*eps)^delta
    edeltat = filter(edelta, filter = c(0, alpha), sides = 1)
    c = omega/(1-sum(beta))
    h = c( h[1:uv], c + filter(edeltat[-(1:uv)], filter = beta,
                               method = "recursive", init = h[uv:1]-c))    
#     if(alpha == 0.1 && flag == TRUE)
#     {
#       print(h)
#       flag = FALSE
#     }
    
    
    
    
    hh = sqrt(abs(h))
    llh = -sum(log(garchDist(z, hh)))
    llh 
  }
  print(garchLLH(params))
  # Step 5: Estimate Parameters and Compute Numerically Hessian:
  fit = nlminb(start = params, objective = garchLLH,
               lower = lowerBounds, upper = upperBounds, control = list(trace=3))
  epsilon = 0.0001 * fit$par
  Hessian = matrix(0, ncol = 4, nrow = 4)
  for (i in 1:4) {
    for (j in 1:4) {
      x1 = x2 = x3 = x4 = fit$par
      x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
      x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
      x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
      x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
      Hessian[i, j] = (garchLLH(x1)-garchLLH(x2)-garchLLH(x3)+garchLLH(x4))/
        (4*epsilon[i]*epsilon[j])
    }
  }
  
  # Step 6: Create and Print Summary Report:
  se.coef = sqrt(diag(solve(Hessian)))
  tval = fit$par/se.coef
  matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(matcoef) = list(names(tval), c(" Estimate",
                                          " Std. Error", " t value", "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  fit$par
}



###############
# Studying the filter process
library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1][1:20]

N = 20; 
# eps = round(rnorm(N), digits = 2) / 10
eps = x - mean(x)
omega = 0.1; alpha = c(0.1); gamma = c(0); beta = c(0.8); delta = 2
u = length(alpha); v = length(beta); uv = max(u,v); h = rep(0.1, uv)

# Conditional Variance filtering - for loop - Wuertz et al. (2006)
for (i in(uv+1):N )
{
    ed = 0
    for (j in 1:u)
    {
        ed = ed+alpha[j]*(abs(eps[i-j])-gamma[j]*eps[i-j])^delta
    }
    h[i] = omega + ed + sum(beta*h[i-(1:v)])
}
  
# Conditional Variance filtering - filter - Wuertz et al. (2006)
edelta = (abs(eps)-gamma*eps)^delta
edeltat = filter(edelta, filter = c(0, alpha), sides = 1)
c = omega/(1-sum(beta))
h = c( h[1:uv], c + filter(edeltat[-(1:uv)], filter = beta,
                             method = "recursive", init = h[uv:1]-c))


# # Conditional Variance filtering - Alternative Approach
# h <- rep(0.1, pq)
# edeltat = 0
# for( i in 1:p) {
#   edelta <- alpha[i]*(abs(z)-gm[i]*z)^delta
#   edeltat = edeltat + edelta[(p-(i-1)):(N-i)]
#   }
# edeltat = c(h[1:p],edeltat)
# c <- omega/(1-sum(beta))
# h <- c( h[1:pq], c + filter(edeltat[-(1:pq)], filter = beta,
#                               > method = "recursive", init = h[q:1]-c))
# hh <- abs(h)^(1/delta)






#################
#################
# Construction of my filter function
#################
# General specifications
  # It must filter only the garch part when needed.
  # It must filter only the arma part when needed.
  # It must filter the both cases when dealing with a combined arma-aparch/garch.
  # It must be a separate function called inside garchFit. 



# Steps on building this function
# First step: Build the filtering function for APARCH models...


filter1.garch11Fit <- function(data, parm)
{
    mu = parm[1]; omega = parm[2]; alpha = parm[3]; beta = parm[4]
    z = (data-mu); Mean = mean(z^2)
    
    # Use Filter Representation:
    e = omega + alpha * c(Mean, z[-length(data)]^2)
    h = filter(e, beta, "r", init = Mean)
    print(h[1])
    hh = sqrt(abs(h))
    cbind(z,hh)  
}




filter.Aparch <- function(
  data,init.Value = NULL,
  p,q, 
  mu, omega, alpha, beta, gamma, delta)
{
  
    # Arguments
    # data: vector with data
    # p,q: model order. aparch(p,q)
    # mu,omega,alpha,beta,gamma,delta: model parameters
    # REMARKS: This function filters, pure 'aparch' models, but 
    # it can deals with other models, such as garch, arch. 
    # All the model parameters must be specified, regardless if they
    # exist. For example:
    # for garch(p,q) with p,q >= 1 make gamma = 0.
    # for arch(p) p >= 1 make q = 1, beta = 0 and gamma = 0.
    
    # Return
    # the series 'z' and 'h'
  
  
    # input 
    # data: vector with data
  
    # error treatment of input parameters
    if( p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
        || length(beta) != q)
        stop("One or more of these conditions were true:
          p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
          || length(beta) != q")
    
    # Initial declaration of variables
    pq = max(p,q)
    z = (data-mu)
    N = length(data)
    Mean.z = mean(abs(z)^delta)
    
    # initializing the time series
    if(is.null(init.Value)) {
        edelta.init <- rep(Mean.z,p)
        h.init <- rep(Mean.z,q)
    } else {
        edelta.init <- rep(init.Value,p)
        h.init <- rep(init.Value,q)  
    }  
    

    edeltat = 0
    for( i in 1:p)
    {
      edelta <- alpha[i]*(c(edelta.init,((abs(z)-gamma[i]*z)^delta)[1:(N-1)]))
      edeltat = edeltat +  edelta[(p-(i-1)):(p+N-i)]
    }
    edeltat = omega + edeltat
    
    h <- filter(edeltat, filter = beta,
                method = "recursive", init = h.init)
    
    if(length(z) != length(h))
      stop("Error in filtering function. length(z) != length(h)")    
    hh <- (abs(h))^(1/delta)
    
    # return
    cbind(z,hh)
}

# Conditional Variance filtering - for loop - Wuertz et al. (2006)
filter.Aparch.Forloop <- function(
    data,h.init = 0.1,
    p,q, 
    mu, omega, alpha, beta, gamma, delta)
{
  
    # Return
    # the series 'z' and 'hh'
    
    
    # input 
    # data: vector with data
    
    # error treatment of input parameters
    if( p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
        || length(beta) != q)
      stop("One or more of these conditions were true:
            p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
            || length(beta) != q")
  
    # Initial declaration of variables
    pq = max(p,q)
    
    z = (data-mu)
    N = length(data)
    Mean.z = mean(abs(z)^delta)
    h = rep(h.init, pq)
    
    # Calculate h[(pq+1):N] recursively
    for (i in (pq+1):N )
    {
        ed = 0
        for (j in 1:p)
        {        
            ed = ed+alpha[j]*(abs(z[i-j])-gamma[j]*z[i-j])^delta
        }
        h[i] = omega + ed + sum(beta*h[i-(1:q)])
    }
    if(length(z) != length(h))
      stop("Error in filtering function. length(z) != length(h)")  
    
    hh <- (abs(h))^(1/delta)
    
    # return
    cbind(z,hh)    
}



################################################################################
# TEST CASES FOF FUNCTIONS:               SPECIFICATION:
#  filter.Aparch                    Test how our filter.Aparh function
#  filter1.garch11Fit               build with matrices operation is related  
# filter.Aparch.Forloop 					  to the other filter functions described
#                                   in Wuertz et al. (2006)
################################################################################

library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1][1:1000]

###############
# Comparison between our function and the filter function from Wurtz.
# garch(1,1)
filter1Result <- filter.Aparch(data = x, p = 1,q = 1,mu = 0.4, omega = 0.1, 
              alpha = c(0.3), beta = c(0.5), gamma = c(0),delta = 2)

filter2Result <- filter1.garch11Fit(x,parm = c(mu = 0.4, omega = 0.1, alpha = 0.3, beta = 0.5))

filter3Result <- filter.Aparch.Forloop(data = x, h.init = 0.5381661^2, p = 1,q = 1,mu = 0.4, omega = 0.1, 
                               alpha = c(0.3), beta = c(0.5), gamma = c(0),delta = 2)

cbind(filter1Result[,1],filter2Result[,1],filter3Result[,1],filter1Result[,2],filter2Result[,2],filter3Result[,2])

# arch(1) == garch(1,0): put p = 1, q = 1 and beta = 0. 
# Notice that the filter.Aparch.Forloop must be equal 
# to the other filter functions after a couple of steps. 
# This is due to the different starting points they have.
filter1Result <- filter.Aparch(data = x, p = 1,q = 1,mu = 0.4, omega = 0.1, 
                               alpha = c(0.3), beta = c(0), gamma = c(0),delta = 2)
filter2Result <- filter1.garch11Fit(x,parm = c(mu = 0.4, omega = 0.1, alpha = 0.3, beta = 0))

filter3Result <- filter.Aparch.Forloop(data = x, h.init = 0.4136527^2, p = 1,q = 1,mu = 0.4, omega = 0.1, 
                                      alpha = c(0.3), beta = c(0), gamma = c(0),delta = 2)
cbind(filter1Result[,1],filter2Result[,1],filter3Result[,1],filter1Result[,2],filter2Result[,2],filter3Result[,2])

# garch(2,2)
filter1Result <- filter.Aparch(data = x, p = 2,q = 2,mu = 0.4, omega = 0.1, 
                               alpha = c(0.3,0.1), beta = c(0.1,0.2), gamma = c(0,0),delta = 2)

filter3Result <- filter.Aparch.Forloop(data = x, h.init = 0.1, p = 2,q = 2,mu = 0.4, omega = 0.1, 
                                       alpha = c(0.3,0.1), beta = c(0.1,0.2), gamma = c(0,0),delta = 2)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])


# arch(15) == garch(15,0)
filter1Result <- filter.Aparch(data = x, p = 15,q = 1,mu = 0.4, omega = 0.1, 
                               alpha = rep(0.03,15), beta = c(0), gamma = rep(0,15),delta = 2)

filter3Result <- filter.Aparch.Forloop(data = x, h.init = 0.1, p = 15,q = 1,mu = 0.4, omega = 0.1, 
                                       alpha = rep(0.03,15), beta = c(0), gamma = rep(0,15),delta = 2)
a <- NULL; 
a <- 0.2
is.null(a)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

# arch(15) == garch(15,0)
filter1Result <- filter.Aparch(data = x, p = 15,q = 1,mu = 0.4, omega = 0.1, 
                               alpha = rep(0.03,15), beta = c(0), gamma = rep(0,15),delta = 2)

filter3Result <- filter.Aparch.Forloop(data = x, h.init = 0.1, p = 15,q = 1,mu = 0.4, omega = 0.1, 
                                       alpha = rep(0.03,15), beta = c(0), gamma = rep(0,15),delta = 2)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

# aparch(4,5)
filter1Result <- filter.Aparch(data = x, p = 4,q = 5,mu = -5, omega = 5, 
                               alpha = rep(0.05,4), beta = c(0.1,0.03,0.02,0.005,0.02), gamma = rep(0.9,4),delta = 1.1)

filter3Result <- filter.Aparch.Forloop(data = x, h.init = 0.1, p = 4,q = 5,mu = -5, omega = 5, 
                                       alpha = rep(0.05,4), beta = c(0.1,0.03,0.02,0.005,0.02), gamma = rep(0.9,4),delta = 1.1)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

# aparch(1,1)
filter1Result <- filter.Aparch(data = x, p = 1,q =1 ,mu = -5, omega = 5, 
                               alpha = rep(0.3), beta = c(0.1), gamma = rep(0),delta = 2)

filter3Result <- filter.Aparch.Forloop(data = x, h.init = 0.1, p = 1,q =1 ,mu = -5, omega = 5, 
                                       alpha = rep(0.3), beta = c(0.1), gamma = rep(0),delta = 2)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

# aparch(2,2)
filter1Result <- filter.Aparch(data = x, p = 2,q =2 ,mu = -5, omega = 5, 
                               alpha = rep(0.1,2), beta = c(0.1,0.3), gamma = rep(0,2),delta = 1.1)

filter3Result <- filter.Aparch.Forloop(data = x, h.init = 0.1, p = 2,q =2 ,mu = -5, omega = 5, 
                                       alpha = rep(0.1,2), beta = c(0.1,0.3), gamma = rep(0,2),delta = 1.1)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

# aparch(1,0)
filter1Result <- filter.Aparch(data = x, p = 1,q =1 ,mu = -5, omega = 5, 
                               alpha = 0.2, beta = 0 , gamma = 0, delta = 20)

filter3Result <- filter.Aparch.Forloop(data = x, h.init = 0.1, p = 1,q =1 ,mu = -5, omega = 5, 
                                       alpha = 0.2, beta = 0 , gamma = 0, delta = 20)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

###############
# Testing different starting points in our filter function

# arch(1) == garch(1,0)
filter1Result <- filter.Aparch(data = x, p = 1,q =1 ,mu = -5, omega = 5, 
                               alpha = 0.2, beta = 0 , gamma = 0, delta = 20)

filter2Result <- filter.Aparch(data = x, p = 1,q =1 ,mu = -5, omega = 5,init.Value = 100, 
                               alpha = 0.2, beta = 0 , gamma = 0, delta = 20)

cbind(filter1Result[,2],filter2Result[,2])
filter1Result[,2]-filter2Result[,2]

# aparch(2,3)
filter1Result <- filter.Aparch(data = x, p = 2,q =3 ,mu = -2, omega = 0.5, 
                               alpha = rep(0.1,2), beta = rep(0.04,3) , gamma = c(0.3,-0.8), delta = 1.2)

filter2Result <- filter.Aparch(data = x, p = 2,q =3 ,mu = -2, omega = 0.5,init.Value = 100, 
                               alpha = rep(0.1,2), beta = rep(0.04,3) , gamma = c(0.3,-0.8), delta = 1.2)

cbind(filter1Result[,2],filter2Result[,2])
filter1Result[,2]-filter2Result[,2]

# garch(59,100): Difficult to achieve the same time series after some steps. We need to pay 
# attention that we are not going to start our time series deliberatelly. We may start it 
# at a place that is representative of the conditional variance.
x <- rnorm(100000)
filter1Result <- filter.Aparch(data = x, p = 59,q =100 ,mu = -2, omega = 0.5, 
                               alpha = rep(5,59), beta = rep(0.04,100) , gamma = rep(0.5,59), delta = 2)

filter2Result <- filter.Aparch(data = x, p = 59,q =100 ,mu = -2, omega = 0.5, init.Value = var((abs(x+2))^2),
                               alpha = rep(5,59), beta = rep(0.04,100) , gamma = rep(0.5,59), delta = 2)

cbind(filter1Result[,2],filter2Result[,2])
filter1Result[,2]-filter2Result[,2]
var((abs(x+2))^2)


