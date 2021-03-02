test_that("unitreg works", {
  set.seed(1212)
  n     <- 200
  betas <- c(1, 2)
  X     <- cbind(1, x1=runif(n))
  eta   <- c(X %*% betas)
  mu    <- exp(eta) / (1 + exp(eta))
  phi   <- exp(6.0)
  y     <- rbeta(n, shape1=mu * phi + 1, shape2=(1-mu) * phi + 1)
  data  <- data.frame(y=y, x1=X[,2])

  m1 <- unitModalReg(y ~ x1, data = data, family = "betaMode")
  m2 <- unitModalReg(y ~ x1, data = data, family = "unitMaxwell")
  m3 <- unitModalReg(y ~ x1, data = data, family = "GPB")
  m4 <- unitModalReg(y ~ x1, data = data, family = "unitGompertz")

  expect_equal(length(m1), 15)
  expect_equal(length(m2), 15)
  expect_equal(length(m3), 15)


})

# solve(
#   numDeriv::hessian(llgpb, x = c(1, 2, log(6.0)), linkinv=make.link("logit")$linkinv, X=X, y=y)
# )








