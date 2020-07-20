#' @title Half-Normal Plot with Simulation Envelopes for \code{unitModalReg} Objects
#'
#' @description Produces a (half-)normal plot from a fitted model object of class \code{\link{unitModalReg}}.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object fitted model object of class \code{\link{unitModalReg}}.
#' @param nsim number of simulations used to compute envelope. Default is 99.
#' @param resid.type type of residuals to be used: The currently options are \code{cox-snell} or \code{quantile}.
#' @param halfnormal logical. If \code{TRUE}, a half-normal plot is produced. If \code{FALSE}, a normal plot is produced.
#' @param plot should the (half-)normal plot be plotted? Default is \code{TRUE}.
#' @param level confidence level of the simulated envelope. Default is 0.95.
#' @param ... additional graphical parameters.
#'
#'
#' @importFrom stats qnorm qexp residuals quantile median rbeta
#' @importFrom graphics plot lines
#'
#' @name hnp
NULL

#' @rdname hnp
#' @export
hnp <- function(object, ...){
  UseMethod("hnp", object)
}


#' @rdname hnp
#' @export
hnp.unitModalReg <- function(object, nsim = 99,
                           halfnormal = TRUE, plot = TRUE, level = 0.95,
                           resid.type = c("cox-snell", "quantile"), ...) {

  if (length(resid.type) > 1) resid.type <- "quantile"

  mu   <- object$fitted.values
  init <- object$coefficients
  p    <- length(init)

  control  <- object$control
  family   <- object$family
  data     <- object$df
  formula  <- object$formula
  link     <- object$link$name
  n        <- nrow(data)
  id       <- as.character(formula)[2]

  if (family != "unitMaxwell") {
    phi <- init[p]
    if (family == "betaMode") phi <- exp(phi)
  }


  ysim <- switch (family,
          "unitMaxwell" = rumaxwell(n = n * nsim, mu = mu),
          "betaMode" = rbeta(n * nsim, shape1 = mu * phi + 1, shape2 = (1-mu) * phi + 1),
          "unitGompertz" = rugompertz(n * nsim, mu = mu, phi = phi),
          "unitGamma" = rugamma(n * nsim, mu = mu, phi = phi),
          "Kumaraswamy" = rkum(n * nsim, mu = mu, phi = phi)
  )

  ysim  <- matrix(ysim, nrow = n, ncol = nsim)

  # Checar se a simulação ta ok!
  # modas <- apply(ysim, 1, function(x) {
  #   d=density(x)
  #   u=d$x
  #   u[which.max(d$y)]
  # } )
  # cbind(modas, mu)[1:10,]


  res_sim <- sapply(1:nsim, function(j) {
    df      <- data
    df[,id] <- ysim[,j]
    obj     <- tryCatch(expr = unitModalReg(
      formula = formula,
      family = family,
      link = link,
      data = df,
      start = init,
      control = control
    ),
    error = function(e) NULL)
    if(!is.null(obj)) {
      residuals(obj, type = resid.type)
    }
    else {
      rep(NaN, n)
    }
  })

  is.na(res_sim) <- sapply(res_sim, is.infinite)

  if (halfnormal) {
    res_obs <- sort(abs(residuals(object, type = resid.type)))
    res_sim <- apply(res_sim, 2, function(x) sort(abs(x), na.last = TRUE))
    if (resid.type == "quantile") {
      res_teo <- qnorm((1:n + n - 1/8) / (2 * n + 0.5))
    }
    if (resid.type == "cox-snell") {
      res_teo <- qexp((1:n + n - 1/8) / (2 * n + 0.5))
    }
  }
  else {
    res_obs <- sort(residuals(object, type = resid.type))
    res_sim <- apply(res_sim, 2, function(x) sort(x, na.last = TRUE))
    if (resid.type == "quantile") {
      res_teo <- qnorm((1:n - 3 / 8) / (n + 1 / 4))
    }
    if (resid.type == "cox-snell") {
      res_teo <- qexp((1:n - 3 / 8) / (n + 1 / 4))
    }
  }

  alpha   <- (1 - level)/2
  res_lwr <- apply(res_sim, 1, quantile, probs = alpha, na.rm = T)
  res_upr <- apply(res_sim, 1, quantile, probs = 1 - alpha, na.rm = T)
  res_mid <- apply(res_sim, 1, median, na.rm = T)

  if (plot) {
    Ry <- c(min(res_lwr), max(res_upr))
    Rx <- range(res_teo)

    plot(
      x = res_teo,
      y = res_obs,
      xlab = 'Theoretical quantiles',
      ylab = 'Residuals',
      xlim = Rx,
      ylim = Ry,
      bty = 'o',
      pch = 3,
      ...
    )
    lines(x = res_teo, y = res_lwr)
    lines(x = res_teo, y = res_upr)
    lines(x = res_teo, y = res_mid, lty = 2)
  }

  list(
    res_obs = res_obs,
    res_teo = res_teo,
    lower   = res_lwr,
    median  = res_mid,
    upper   = res_upr
  )
}
