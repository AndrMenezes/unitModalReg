test_that("summary works", {
  set.seed(1212)
  n     <- 200
  betas <- c(1, 2)
  X     <- cbind(1, x1=runif(n))
  eta   <- c(X %*% betas)
  mu    <- exp(eta) / (1 + exp(eta))
  phi   <- exp(6.0)
  y     <- rbeta(n, shape1=mu * phi + 1, shape2=(1-mu) * phi + 1)
  data  <- data.frame(y=y, x1=X[,2])

  m <- unitModalReg(y ~ x1, data = data, family = "betaMode")
  s <- summary(m)
  expect_equal(names(s), c("coeftable", "loglik", "correlation", "call",
                           "link", "family", "data"))

})
