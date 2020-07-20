test_that("residuals works", {

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
  m3 <- unitModalReg(y ~ x1, data = data, family = "unitGompertz")


  qa1 <- residuals(m1, type = "quantile")
  qa2 <- residuals(m2, type = "quantile")
  qa3 <- residuals(m3, type = "quantile")

  expect_equal(length(qa1), n)
  expect_equal(length(qa2), n)
  expect_equal(length(qa3), n)

  cs1 <- residuals(m1, type = "cox-snell")
  cs2 <- residuals(m2, type = "cox-snell")
  cs3 <- residuals(m3, type = "cox-snell")

  expect_equal(length(cs1), n)
  expect_equal(length(cs2), n)
  expect_equal(length(cs3), n)



})
