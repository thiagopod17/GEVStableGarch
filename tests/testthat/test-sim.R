test_that("GARCH(1,1) with GAT do not give NaN", {
  
  skip_on_cran()
  
  set.seed(123)
  
  true <- c(mu = 0, omega = 0.4, alpha1 = 0.2, beta1 = 0.2,
            skew = 0.5, shape1 = 2, shape2 = 2)
  
  
  spec.gat = GEVStableGarch::gsSpec(model = list(mu = true["mu"], omega = true["omega"],
                                                 alpha = true["alpha1"], beta = true["beta1"], 
                                                 skew = true["skew"],
                                                 shape = c(true["shape1"], true["shape2"])), 
                                    cond.dist = "gat")
  sim.last.value = as.numeric(GEVStableGarch::gsSim(spec = spec.gat, n = 10000)[,1])[10000]
  expect_false(is.nan(sim.last.value))
  
})







