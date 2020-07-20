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
#' @importFrom stats pbeta
#'

#' @rdname residuals.unitModalReg
#' @export
residuals.unitModalReg <- function(object, type = c("cox-snell", "quantile"), ...) {

  if (length(type) > 1) type = "quantile"


  mu  <- object$fitted.values
  y   <- object$data$y
  phi <- object$coefficients
  phi <- phi[length(phi)]

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

  res
}
