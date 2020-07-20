test_that("hnp works", {

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

  qa1 <- hnp(m1, plot = TRUE, halfnormal = FALSE)
  qa2 <- hnp(m2, plot = TRUE, halfnormal = FALSE)
  qa3 <- hnp(m3, plot = TRUE, halfnormal = FALSE)

  expect_equal(length(qa1), 5)
  expect_equal(length(qa2), 5)
  expect_equal(length(qa3), 5)



})
