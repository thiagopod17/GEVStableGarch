test_that("Passing invalid input parameters for dgat", {
  # sd cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = -1, nu = 2, d = 3, xi = 1))
  
  # d cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 0.1, d = -1, xi = 1))
  
  # xi cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 1, d = 0.1, xi = 0))
  
  # nu cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 0, d = 0.1, xi = 1))
})


test_that("Passing invalid input parameters for dgat", {
  # sd cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = -1, nu = 2, d = 3, xi = 1))
  
  # d cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 0.1, d = -1, xi = 1))
  
  # xi cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 1, d = 0.1, xi = 0))
  
  # nu cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 0, d = 0.1, xi = 1))
})


test_that("Integration of dgat is approximately 1", {
  expect_equal(integrate(dgat , lower = -Inf, upper = Inf, 
                         mean = 0.3, sd = 1, nu = 1,
                         d = 2, xi = 1)$value, 1, tolerance = 1e-5)
})


test_that("Distribution function behaviour close to -Inf and Inf", {
  expect_equal(pgat(-1000, mean = 0.3, sd = 1, nu = 1,
                    d = 2, xi = 1), 0, tolerance = 1e-5)
  
  expect_equal(pgat(1000, mean = 0.3, sd = 1, nu = 1,
                    d = 2, xi = 1), 1, tolerance = 1e-5)
})


test_that("Distribution and Density property: integral( pgat, a, b) = pgat(b) - pgat(a)", {
  a = -3
  b = 5
  prob_with_integral = integrate(dgat , lower = a, upper = b, 
                                 mean = 0.3, sd = 2, nu = 1,
                                 d = 2, xi = 1)$value
  pgat_at_b = pgat(b, mean = 0.3, sd = 2, nu = 1,
                   d = 2, xi = 1)
  pgat_at_a = pgat(a, mean = 0.3, sd = 2, nu = 1,
                   d = 2, xi = 1)
  
  expect_equal(prob_with_integral, pgat_at_b - pgat_at_a, tolerance = 1e-5)
})


test_that("Distribution and Density property: integral( pgat, a, b) = pgat(b) - pgat(a)", {
  
  prob_value = 0.344
  quantile_value = qgat(prob_value, mean = 0.3, sd = 2, nu = 1,
       d = 2, xi = 1)
  prob_value_with_distribution = pgat(quantile_value, mean = 0.3, sd = 2, nu = 1,
       d = 2, xi = 1)
  expect_equal(prob_value, prob_value_with_distribution, tolerance = 1e-5)

})



test_that("Expectation via Monte Carlo matches numerical integration", {
  
  set.seed(123)
  
  # Parameters
  a  <- -Inf
  b  <- Inf
  mean <- 0.3
  sd   <- 2
  nu   <- 1
  d    <- 2
  xi   <- 1
  
  # Numerical integration for E[X]
  integral_mean <- integrate(
    function(x) {
      x * dgat(
        x,
        mean = mean,
        sd   = sd,
        nu   = nu,
        d    = d,
        xi   = xi
      )
    },
    lower = a,
    upper = b
  )$value
  
  # Monte Carlo approximation
  n <- 1e5
  samples <- rgat(
    n,
    mean = mean,
    sd   = sd,
    nu   = nu,
    d    = d,
    xi   = xi
  )
  
  mc_mean <- mean(samples)
  
  # Expect close agreement
  expect_equal(
    mc_mean,
    integral_mean,
    tolerance = 0.05
  )
})



