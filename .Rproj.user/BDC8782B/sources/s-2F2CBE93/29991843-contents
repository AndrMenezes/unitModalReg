#' @title Control Parameters for unit modal regression models
#'
#' @description Parameters that control fitting of unit modal regression models using in \code{unitModalReg}.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param method characters string specifying the method argument passed to \code{optim}.
#' @param maxit integer specifying the \code{maxit} argument (maximal number of iterations) passed to \code{optim}.
#' @param reltol relative convergence tolerance passed to \code{optim}.
#' @param ... arguments passed to \code{optim}.
#'

#' @rdname unitModalReg.control
#' @export
unitModalReg.control <- function(method = "BFGS", maxit = 5000, reltol = 1e-8, ...) {
  out <- list(
    method = method,
    maxit = maxit,
    reltol = reltol
  )
  out <- c(out, list(...))
  if (!is.null(out$fnscale)) warning("fnscale must not be modified")
  out$fnscale <- 1
  out
}

#' @title Unit Modal Regression Models
#'
#' @description Fit several unit modal regression models.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param formula symbolic description of the quantile model.
#' @param data data.frame contain the variables in the model.
#' @param family specify the distribution family name.
#' @param link character specification of the link function in the modal regression. Default is \code{logit}.
#' @param start an optional vector with starting values for all parameters
#' @param control a list of control arguments specified via \code{unitModalReg.control}.
#' @param y numeric vector of response variable.
#' @param X numeric matrix. Design matrix for the modal regression.
#'
#'
#' @importFrom stats optim make.link model.frame model.matrix model.response
#' @importFrom Formula as.Formula
#' @name unitModalReg
NULL

#' @rdname unitModalReg
#' @export
unitModalReg <- function(formula, data, family, link, start = NULL, control = unitModalReg.control()) {


  if (missing(link)) link <- "logit"

  # building model matrices
  ff <- Formula::as.Formula(formula)
  mf <- model.frame(ff, data = data)
  y  <- model.response(mf, "numeric")
  X  <- model.matrix(ff, data = data, rhs = 1)
  p  <- ncol(X)

  ## check response variable
  if(!(min(y) > 0 & max(y) < 1)) {
    stop("invalid dependent variable, all observations must be in (0, 1)")
  }

  fit <- unitModalReg.fit(
    y = y,
    X = X,
    family = family,
    link = link,
    start = start,
    control = control
  )

  if (!fit$converged) {
    stop("optimization failed to converge!")
  }

  vcov <- tryCatch(solve(fit$hessian), error = function(e) NULL)
  if (is.null(vcov)) {
    llike <- switch (family,
                  "unitMaxwell" = llumaxwell,
                  "betaMode" = llbetamode,
                  "unitGompertz" = llugompertz,
                  "unitGamma" = llugamma,
                  "Kumaraswamy" = llkum,
                  "GPB" = llgpb
    )
    numDeriv::hessian(llike, x = fit$par, linkinv=make.link(link)$linkinv, X=X, y=y)
  }

  # fitted values
  parms         <- fit$par[1:p]
  linkobj.mu    <- make.link(link)
  fitted.values <- linkobj.mu$linkinv(X %*% parms)

  # Output
  out <- list(
    call             = match.call(),
    formula          = ff,
    control          = control,
    link             = linkobj.mu,
    family            = family,
    loglik           = -fit$value,
    vcov             = vcov,
    coefficients     = fit$par,
    fitted.values    = c(fitted.values),
    nobs             = length(y),
    npar             = length(fit$par),
    df.residual      = length(y) - length(fit$par),
    data             = list(X = X, y = y),
    df               = data
  )
  class(out) <- "unitModalReg"
  out
}

#' @rdname unitModalReg
#' @export
unitModalReg.fit <- function(y, X, family, link, start, control = unitModalReg.control()) {

  llike <- switch (family,
    "unitMaxwell" = llumaxwell,
    "betaMode" = llbetamode,
    "unitGompertz" = llugompertz,
    "unitGamma" = llugamma,
    "Kumaraswamy" = llkum,
    "GPB" = llgpb
  )

  n <- length(y)
  p <- ncol(X)

  linkobj <- stats::make.link(link)
  method  <- control$method
  control$method <- NULL

  # Initial guess
  if (is.null(start)) {
    # For beta
    ystar        <- linkobj$linkfun(y)
    reg_ini      <- suppressWarnings(stats::lm.fit(x = X, y = ystar))
    start        <- reg_ini$coefficients
    names(start) <- colnames(X)

    # For phi/gamma
    if (family != "unitMaxwell") {
      yhat  <- linkobj$linkinv(reg_ini$fitted.values)
      p     <- length(reg_ini$coefficients)
      phi   <- est_phi(family, y, yhat, p)
      start <- c(start, phi = phi)
    }
  }

  # Maximization
  opt <- optim(
    par     = start,
    fn      = llike,
    method  = method,
    hessian = TRUE,
    control = control,
    X       = X,
    y       = y,
    linkinv = linkobj$linkinv
  )
  ## check if the optim converged
  if (opt$convergence > 0) {
    opt$converged <- FALSE
    warning("optimization failed to converge")
  } else {
    opt$converged <- TRUE
  }
  opt
}
