#' @title Residuals Method for \code{unitModalReg} Objects
#'
#' @description Extract various types of residuals from unit regression models: cox-snell and quantile residuals for now.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object fitted model object of class \code{\link{unitModalReg}}.
#' @param type character indicating type of residuals.
#' @param ... currently not used.
#'
#'
#' @importFrom stats pbeta predict
#' @importFrom numDeriv grad
#'

#' @rdname residuals.unitModalReg
#' @export
residuals.unitModalReg <- function(object, type = c("cox-snell", "quantile",
                                                    "pearson", "working",
                                                    "response"), ...) {

  if (length(type) > 1) type = "quantile"


  mu <- object$fitted.values
  y <- object$data$y
  parms <- object$coefficients
  phi <- parms[length(parms)]

  if (type %in% c("quantile", "cox-snell")) {
    Fy <- switch (object$family,
                  "unitMaxwell" = pumaxwell(y, mu=mu),
                  "betaMode" = pbeta(y, shape1=mu * phi + 1, shape2=(1-mu) * phi + 1),
                  "unitGompertz" = pugompertz(y, mu=mu, phi=phi),
                  "unitGamma" = pugamma(y, mu=mu, phi=phi),
                  "Kumaraswamy" = pkum(y, mu=mu, phi=phi)
    )
    res <- switch (type,
                   "cox-snell" = -log(1-Fy),
                   "quantile" =  qnorm(Fy)
    )
  }
  if (type == "response") {
    res <- y - mu
  }
  if (type %in% c("working", "partial")) {
    link <- object$link
    eta <- link$linkfun
    d_eta_mu <- numDeriv::grad(func = eta, x = mu)
    eta_hat <- eta(mu)
    r <- (y - mu) * d_eta_mu

    res <- switch (type,
      "working" = r,
      "partial" = r + predict(object, type = "terms")
    )

  }
  if (type == "pearson") {
    Vary <- compute_var(family = object$family, mu = mu, phi = phi)
    res <- (y - mu) / sqrt(Vary)
  }

  return(res)
}
