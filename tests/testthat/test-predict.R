test_that("predict works", {
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
  pred_response_1 <- predict(m, type = "response")
  pred_response_2 <- predict(m, type = "response", interval = TRUE)
  pred_response_3 <- predict(m, type = "response", interval = TRUE,
                             se.fit = TRUE)

  expect_equal(ncol(pred_response_1), 1, label = "predict type response")
  expect_equal(ncol(pred_response_2), 3,
               label = "predict type response with interval")
  expect_equal(ncol(pred_response_3), 4,
               label = "predict type response with interval and se")

  pred_terms <- predict(m, type = "terms")
  expect_equal(ncol(pred_terms), 1,
               label = "predict type terms with one covariate")

  data$x2 <- rbinom(n, size = 1, prob = 0.5)
  m <- unitModalReg(y ~ x1 + x2, data = data, family = "betaMode")
  pred_terms <- predict(m, type = "terms")
  expect_equal(ncol(pred_terms), 2,
               label = "predict type terms with two covariate")

  # Check predict with newdata
  data2predict <- data.frame(x1 = runif(4), x2 = rexp(4))
  pred_response_1 <- predict(m, newdata = data2predict, type = "response")
  pred_response_2 <- predict(m, newdata = data2predict, type = "response",
                             interval = TRUE)
  pred_response_3 <- predict(m, newdata = data2predict, type = "response",
                             interval = TRUE, se.fit = TRUE)

  expect_equal(ncol(pred_response_1), 3, label = "predict with newdata")
  expect_equal(ncol(pred_response_2), 5,
               label = "predict with newdata and interval")
  expect_equal(ncol(pred_response_3), 6,
               label = "predict with newdata, interval and se")

})
