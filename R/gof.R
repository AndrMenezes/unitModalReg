#' @title Goodness-of-Fit for models of class \code{unitModalReg}
#'
#' @description Compute goodness-of-fit measures and Voung test for models of class \code{unitModalReg}.
#'
#' @author André F. B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param ... pass the models separate by commma to make the comparisons or use the argument \code{lt}.
#' @param lt list of models.
#' @param object1,object2 obected of class \code{unitModalReg} contain the fitted models.
#'
#'
#' @importFrom stats AIC BIC var
#' @importFrom utils combn
#' @name Goodness-of-Fit
NULL

#' @rdname Goodness-of-Fit
#' @export

gof <- function(..., lt=NULL) {

  if (is.null(lt)) lt <- list(...)

  n <- length(lt[[1]]$data$y)

  info_crit <- sapply(lt, function(x) {
    aic <- AIC(x)
    bic <- BIC(x)
    hqic <- AIC(x, k = 2 * log(log(n)))
    c(aic=aic,bic=bic,hqic=hqic)
  })
  colnames(info_crit) <- sapply(lt, get, x="family")

  objs <- combn(lt, 2,simplify=FALSE)

  voung <- lapply(1:length(objs), function(j) {
    out <- matrix(unlist(voung.test(objs[[j]][[1]], objs[[j]][[2]])[c("statistic", "p.value")]), ncol=1)
    rownames(out) <- c("statistic", "p.value")
    colnames(out) <- paste0(objs[[j]][[1]]$family, " versus ", objs[[j]][[2]]$family)
    out
  })
  voung <- do.call("cbind", voung)

  out <- list(
    call=match.call(),
    voung=voung,
    info_crit=info_crit
  )
  class(out) <- "gof"
  out
}

#' @rdname Goodness-of-Fit
#' @export
voung.test <- function(object1, object2) {

  y <- object1$data$y
  X <- object1$data$X
  n <- length(y)

  ll_1 <- switch(
    object1$family,
    "unitMaxwell" = llumaxwell,
    "betaMode" = llbetamode,
    "GPB" = llgpb,
    "unitGompertz" = llugompertz,
    "unitGamma" = llugamma,
    "Kumaraswamy" = llkum
  )
  ll_2 <- switch(
    object2$family,
    "unitMaxwell" = llumaxwell,
    "betaMode" = llbetamode,
    "GPB" = llgpb,
    "unitGompertz" = llgpb,
    "unitGamma" = llugamma,
    "Kumaraswamy" = llkum
  )
  par_1 <- object1$coefficients
  par_2 <- object2$coefficients
  linkinv_1 <- object1$link$linkinv
  linkinv_2 <- object2$link$linkinv

  ll_1 <- ll_1(par_1, linkinv_1, X, y, sum=FALSE)
  ll_2 <- ll_2(par_2, linkinv_1, X, y, sum=FALSE)

  om2     <- (n - 1) / n * var(ll_1 - ll_2)
  lr      <- sum(ll_1 - ll_2)
  Tstat   <- (1/sqrt(n)) * lr/sqrt(om2)
  pvalue  <- pnorm(q = abs(Tstat), lower.tail = F)
  names(Tstat) <- "T_LR"


  out <- list(
    statistic = Tstat,
    p.value = pvalue,
    method = "Vuong likelihood ratio test for non-nested models",
    data = paste0(object1$family, " versus ", object2$family)
  )

  class(out) <- "htest"
  out
}
