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

  expect_equal(length(qa1), n,
               label = "quantile residuals for beta modal regression")
  expect_equal(length(qa2), n,
               label ="quantile residuals for unitMaxwell modal regression")
  expect_equal(length(qa3), n,
               label ="quantile residuals for unitGompertz modal regression")

  cs1 <- residuals(m1, type = "cox-snell")
  cs2 <- residuals(m2, type = "cox-snell")
  cs3 <- residuals(m3, type = "cox-snell")

  expect_equal(length(cs1), n,
               label = "cox-snell residuals for beta modal regression")
  expect_equal(length(cs2), n,
               label = "cox-snell residuals for unitMaxwell modal regression")
  expect_equal(length(cs3), n,
               label = "cox-snell residuals for unitGompertz modal regression")


  working <- residuals(m1, type = "working")
  expect_equal(length(working), n,
               label = "working residuals for beta modal regression")

  response <- residuals(m1, type = "response")
  expect_equal(length(response), n,
               label = "response residuals for beta modal regression")

  pearson <- residuals(m1, type = "pearson")
  expect_equal(length(pearson), n,
               label = "pearson residuals for beta modal regression")

  partial <- residuals(m1, type = "partial")
  expect_equal(colnames(partial), "x1",
               label = "partial residuals with one covariate")

  # Adding one more covariate
  data$x2 <- rexp(n)
  m <- unitModalReg(y ~ x1 + x2, data = data, family = "betaMode")
  partial <- residuals(m, type = "partial")
  expect_equal(colnames(partial), c("x1", "x2"),
               label = "partial residuals with one covariate")


})
