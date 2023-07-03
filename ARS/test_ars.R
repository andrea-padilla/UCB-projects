set.seed(8)
n <- 10000

test_that("Correct sampling for standard normal distribution", {

  # use Z = (X_1_bar - X_2_bar) / sqrt(sd_X1^2 + sd_X2^2)
  m = mean(ars(dnorm, 10000))
  s_d = sd(ars(dnorm, 10000))
  # the mean of the standard normal distribution is 0
  # the sd of the standard normal distribution is 1
  Z = (m - 0) / sqrt(1 + s_d^2)
  expect_lte(abs(Z), 1.96)

})

test_that("Correct sampling for unnormalized standard normal distribution", {
  unnorm_gaus <- function(x) {2*exp(-1/2*x^2)/sqrt(2*pi)}

  # use Z = (X_1_bar - X_2_bar) / sqrt(sd_X1^2 + sd_X2^2)
  m = mean(ars(unnorm_gaus, 10000))
  s_d = sd(ars(unnorm_gaus, 10000))
  # the mean of the standard normal distribution is 0
  # the sd of the standard normal distribution is 1
  Z = (m - 0) / sqrt(1 + s_d^2)
  expect_lte(abs(Z), 1.96)

})


test_that("Correct sampling for beta distribution", {
  beta_func <- function(x) {dbeta(x, shape1 = 3, shape2 = 4)} #Beta(3, 4)

  # use Z = (X_1_bar - X_2_bar) / sqrt(sd_X1^2 + sd_X2^2)
  m = mean(ars(beta_func, 10000, absc = c(0.2, 0.6), lower = 0, upper = 1))
  s_d = sd(ars(beta_func, 10000, absc = c(0.2, 0.6), lower = 0, upper = 1))
  # the mean of the distribution Beta(3, 4) is 3/7
  # the var of the distribution Beta(3, 4) is 3/98
  Z = (m - 3/7) / sqrt(3/98 + s_d^2)
  expect_lte(abs(Z), 1.96)

})

test_that("Correct sampling for Gamma distribution", {
  gamma_func <- function(x) {dgamma(x, shape = 3, rate = 4)} # Gamma(3, 4)

  m <- mean(ars(gamma_func, 10000, absc = c(0.1, 5), lower = -Inf, upper = Inf))
  s_d <- sd(ars(gamma_func, 10000, absc = c(0.1, 5), lower = -Inf, upper = Inf))
  # the mean of the distribution Gamma(3, 4) is 3/4
  # the var of the distribution Gamma(3, 4) is 3/16
  Z = (m - 3/4) / sqrt(3/16 + s_d^2)
  expect_lte(abs(Z), 1.96)

})


test_that("Correct sampling for logistic distribution", {
  logis_func <- function(x) {dlogis(x, location = 1, scale = 2.5)} #Logis(1, 2.5)

  # use Z = (X_1_bar - X_2_bar) / sqrt(sd_X1^2 + sd_X2^2)
  m = mean(ars(logis_func, 10000, absc = c(-2, 3), lower = -Inf, upper = Inf))
  s_d = sd(ars(logis_func, 10000, absc = c(-2, 3), lower = -Inf, upper = Inf))
  # the mean of the distribution Logis(1, 2.5) is 1
  # the var of the distribution Logis(1, 2.5) is 2.5^2*pi^2 / 3
  Z = (m - 1) / sqrt(2.5^2 * pi^2 / 3 + s_d^2)
  expect_lte(abs(Z), 1.96)

})


test_that("Correct sampling for laplace distribution", {
  laplace_func <- function(x) {1/(2*1.5) * exp(-abs(x-3)/1.5)} #Laplace(3, 1.5)

  # use Z = (X_1_bar - X_2_bar) / sqrt(sd_X1^2 + sd_X2^2)
  m = mean(ars(laplace_func, 10000, absc = c(1, 4), lower = -5, upper = 5))
  s_d = sd(ars(laplace_func, 10000, absc = c(1, 4), lower = -5, upper = 5))
  # the mean of Laplace(3, 1.5) is 3
  # the var of Laplace(3, 1.5) is 2 * 1.5^2
  Z = (m - 3) / sqrt(2 * 1.5^2 + s_d^2)
  expect_lte(abs(Z), 1.96)

})

test_that("Correct sampling for exponential distribution", {
  exp_adjust <- function(x) {dexp(x, rate = 2.5)} #Just changing to be rate = 2.5

  # use Z = (X_1_bar - X_2_bar) / sqrt(sd_X1^2 + sd_X2^2)
  m = mean(ars(exp_adjust, n, absc = c(1, 4), lower = 0, upper = Inf))
  s_d = sd(ars(exp_adjust, n, absc = c(1, 4), lower = 0, upper = Inf))
  # the mean of the exponential distribution with rate 2.5 is 1/2.5
  # the var of the exponential distribution with rate 2.5 is 1/(2.5^2)
  Z = (m - 1/2.5) / sqrt(1/(2.5^2) + s_d^2)
  expect_lte(abs(Z), 1.96)

})

test_that("Correct sampling for uniform distribution", {
  unif_func <- function(x) {dunif(x, min = 2, max = 5)} #Unif(2,5)

  # use Z = (X_1_bar - X_2_bar) / sqrt(sd_X1^2 + sd_X2^2)
  m = mean(ars(unif_func, n, absc = c(2.5, 4), lower = 2, upper = 5))
  s_d = sd(ars(unif_func, n, absc = c(2.5, 4), lower = 2, upper = 5))
  # the mean of the uniform(2, 5) distribution is (2+5)/2
  # the var of the uniform(2, 5) distribution is (5-2)^2/12
  Z = (m - (2+5)/2) / sqrt((5-2)^2/12 + s_d^2)
  expect_lte(abs(Z), 1.96)

})

### Log-Concave Checks
# Cauchy(2, 4)
cauchy_func <- function(x) {dcauchy(x, location = 2, scale = 4)}
expect_error(ars(cauchy_func, n, absc = c(-3, 3), lower = -Inf, upper = Inf)) #This should fail, Cauchy is non-log-concave.

# t(5)
t_func <- function(x) {dt(x, df = 5)}
expect_error(ars(t_func, n, absc = c(-3, 3), lower = -Inf, upper = Inf)) #This should fail, t-distribution is non-log-concave.

### Correct Abscissae Points Checks
# Normal(0,1) with Abscissae = c(-3, -2). Both initial abscissae points have positive slope, so it fails because unbounded from above.
expect_error(ars(dnorm, n, absc = c(-3, -2))) #unbounded from both sides
expect_error(ars(dnorm, n, absc = c(-3, -2), lower = -10)) #unbounded from above
expect_error(ars(dnorm, n, absc = c(-2, -1), upper = 10), NA) #unbounded from below, not an issue

# Logistic(1, 2) with Abscissae = c(2, 3). Both initial abscissae points have negative slope, so it fails because unbounded from below.
expect_error(ars(logis_func, n, absc = c(2, 3)))
expect_error(ars(logis_func, n, absc = c(2, 3), upper = 10))
expect_error(ars(logis_func, n, absc = c(2, 3), lower = -10), NA) #unbounded from above, not an issue

