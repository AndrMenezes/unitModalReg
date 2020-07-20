#' @title Methods for \code{unitModalReg} Objects
#'
#' @description Methods for extracting information from fitted models objects of class \code{\link{unitModalReg}}.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param object,x fitted model object of class \code{\link{unitModalReg}}.
#' @param digits  minimal number of _significant_ digits
#' @param correlation logical; if \code{TRUE}, the correlation matrix of
#'   the estimated parameters is returned and printed. Default is \code{FALSE}.
#' @param parm a specification of which parameters are to be given confidence intervals,
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param ... additional argument(s) for methods. Currently not used.
#'
#' @importFrom stats qnorm pnorm cov2cor coef vcov printCoefmat
#'
#' @name methods-unitModalReg
NULL


#' @rdname methods-unitModalReg
#' @export

print.unitModalReg <- function(x, digits = 4, ...)
{
  p <- ncol(x$data$X)

  cat("\n", paste0(x$family, " regression model"), sep = "")

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Mu coefficients (mode model with ", x$link$name, " link): \n", sep = "")

  print.default(FF(x$coefficients[1:p], Digits = digits), print.gap = 2, quote = FALSE)

  cat("\n")

  if (x$family != "unitMaxwell") {
    cat("Model with constant shape parameter:", "\n", sep = "")

    print.default(format(x$coefficients[-(1:p)], digits = digits), print.gap = 2, quote = FALSE)

    cat("\n")
  }
  invisible(x)
}


# Summary -----------------------------------------------------------------

#' @rdname methods-unitModalReg
#' @export
summary.unitModalReg <- function(object, correlation = FALSE, ...) {

  estimates <- object$coefficients
  stderror  <- sqrt(diag(object$vcov))
  zvalue    <- estimates/stderror
  pvalue    <- 2 * pnorm(-abs(zvalue))
  table     <- cbind("Estimate"    = estimates,
                     "Std. Error"  = stderror,
                     "Z value"     = zvalue,
                     "Pr(>|z|)"    = pvalue)
  if (correlation) {
    correlation <- cov2cor(object$vcov)
  }

  out <- list(coeftable   = table,
              loglik      = object$loglik,
              correlation = correlation,
              call        = object$call,
              tau         = object$tau,
              link.mu     = object$link$name,
              link.phi    = object$link.phi$name,
              model       = object$model,
              data        = object$data)
  class(out) <- "summary.unitModalReg"
  out
}


# Print output summary ----------------------------------------------------

#' @export
print.summary.unitModalReg <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  p <- ncol(x$data$X)

  cat("\n Wald-tests for ", x$family, " regression model", "\n" ,sep = "")

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Mu coefficients (mode model with ", x$link$name, " link): \n", sep = "")
  printCoefmat(x$coeftable[1:p, , drop = FALSE], digits = digits, has.Pvalue = TRUE)
  cat("\n")

  if (x$family != "unitMaxwell") {
    cat("Model with constant shape:", "\n", sep = "")
    printCoefmat(x$coeftable[-(1:p), , drop = FALSE], digits = digits, has.Pvalue = TRUE)
    cat("\n")
  }

  if (is.matrix(x$correlation)) {
    cat("Correlation of coefficients:", "\n", sep = "")
    corr <- x$correlation
    corr <- format(round(corr, 2L), nsmall = 2L, digits = digits)
    corr[!lower.tri(corr)] <- ""
    print(corr[-1, -ncol(corr), drop = FALSE], quote = FALSE)
    cat("\n")
  }

  invisible(x)
}

# coef function -----------------------------------------------------------

#' @rdname methods-unitModalReg
#' @export
coef.unitModalReg <- function(object, ...) {
  if (!missing(...)) {
    warning("Extra arguments discarded")
  }
  object$coefficients
}


# vcov function -----------------------------------------------------------

#' @rdname methods-unitModalReg
#' @export
vcov.unitModalReg <- function(object, ...) {
  if (!missing(...)) {
    warning("Extra arguments discarded")
  }
  object$vcov
}


# logLik function ---------------------------------------------------------

#' @rdname methods-unitModalReg
#' @export
logLik.unitModalReg <- function(object, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  ll <- object$loglik
  attr(ll, "df")   <- object$npar
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  ll
}

# confint function --------------------------------------------------------

#' @rdname methods-unitModalReg
#' @export
confint.unitModalReg <- function(object, parm, level = 0.95, ...)
{
  cf <- coef(object)
  ses <- sqrt(diag(vcov(object)))
  pnames <- names(ses)
  if (missing(parm)) {
    parm <- pnames
  }
  else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- format.perc(a, 3)
  ci <- array(NA_real_, dim = c(length(parm), 2L),
              dimnames = list(parm, pct))
  ci[] <- cf[parm] + ses[parm] %o% fac
  ci
}

#' @rdname methods-unitModalReg
#' @export
fitted.unitModalReg <- function(object, ...) {

  if (!missing(...)) {
    warning("Extra arguments discarded")
  }
  object$fitted.values
}

# Print output of GoF Test ----------------------------------------------

#' @export
print.gof <- function(x, ...) {

  vo <- FF(t(x$voung), 3)
  cr <- FF(x$info_crit, 3)

  cat("\nGoodness-of-Fit for modal regression models \n", sep = "")

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Voung tests for all model combinations: \n")

  print(vo, quote = F)

  cat("\n\n")

  cat("Information criterions of modal regression models: \n")
  print(cr, quote = F)

}



