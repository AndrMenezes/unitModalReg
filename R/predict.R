#' @title Prediction Method for \code{unitModalReg} Objects
#'
#' @description Extract various types of predictions from fitted models objects
#' of class \code{\link{unitModalReg}}.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object fitted model object of class \code{\link{unitModalReg}}.
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the original observations are used.
#' @param type the type of prediction required. The default is on the scale of
#' the response variable. The "\code{terms}" option returns a matrix giving
#' the fitted values of each term in the model formula on the linear predictor
#' scale.
#' @param interval type of interval desired. The options are \code{none} and \code{confidence}.
#' @param level coverage probability for the confidence intervals. Default is \code{0.95}.
#' @param se.fit logical. If \code{TRUE} return the asymptotic standard errors.
#' @param ... currently not used.
#'
#'
#' @importFrom stats qnorm terms formula
#' @importFrom Formula as.Formula


#' @rdname predict.unitModalReg
#' @export

predict.unitModalReg <- function(object, newdata = NULL, type,
                                 interval = FALSE, level = 0.95,
                                 se.fit = FALSE, ...) {
  parms <- object$coefficients
  betas <- parms[-length(parms)]
  p <- length(betas)
  phi <- parms[length(parms)]
  vcov <- object$vcov
  ff <- Formula::as.Formula(object$formula)
  if (missing(type))
    type <- "response"
  if (is.null(newdata)) {
    X <- object$data$X
  }
  else {
    X <- model.matrix(ff, newdata, rhs = 1)
  }
  if (type == "response") {
    linkobj <- object$link
    vcov.beta <- vcov[1:p, 1:p]
    g.mu.hat <- as.numeric(X %*% betas)
    mu.hat <- linkobj$linkinv(g.mu.hat)
    out <- cbind(fit = mu.hat)
    if (interval) {
      J <- linkobj$mu.eta(mu.hat) * X
      variance <- tcrossprod(J %*% vcov.beta, J)
      stderror <- sqrt(diag(variance))
      qa <- qnorm(1 - level/2)
      out <- cbind(out, lower = out[, "fit"] - qa *
                     stderror, upper = out[, "fit"] + qa * stderror)
    }
    if (se.fit) {
      out <- cbind(out, se.fit = stderror)
    }
  }
  else if (type == "terms") {
    X1 <- sweep(X, 2L, colMeans(X))
    terms_betas <- labels(terms(ff))
    out <- t(betas * t(X1))[, -1]
    out <- as.data.frame(out)
    names(out) <- terms_betas
  }
  if (is.null(newdata)) {
    return(out)
  }
  else {
    return(as.data.frame(cbind(newdata, out)))
  }
}
