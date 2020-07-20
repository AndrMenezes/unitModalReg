#' @title Compute variance
#'
#' @description  Compute analitycal or Monte Carlo variance for unit Modal distributions
#'
#' @author André Felipe B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#' @param family character specifying the model. Currently are: \code{betaMode}, \code{unitGompertz},
#' \code{unitGamma} and \code{Kumaraswamy}
#' @param mu mode parameter
#' @param phi nuance parameter
#' @param MC logical indicating if use Monte Carlo simulations to compute variance.
#'
#'
#' @examples
#' mu <- 0.5
#' phi <- 3.0
#'
#' V_beta   <- compute_var("betaMode", mu, phi)
#' V_ugtz   <- compute_var("unitGompertz", mu, phi)
#' V_ugamma <- compute_var("unitGamma", mu, phi)
#' V_kum    <- compute_var("Kumaraswamy", mu, phi)
#'
#' cbind(mu=mu, phi=phi, V_beta, V_ugtz, V_ugamma, V_kum)
#'
#' V_beta   <- compute_var("betaMode", mu, phi, MC=TRUE)
#' V_ugtz   <- compute_var("unitGompertz", mu, phi, MC=TRUE)
#' V_ugamma <- compute_var("unitGamma", mu, phi, MC=TRUE)
#' V_kum    <- compute_var("Kumaraswamy", mu, phi, MC=TRUE)
#'
#' cbind(mu=mu, phi=phi, V_beta, V_ugtz, V_ugamma, V_kum)
#'

#' @export
compute_var <- function(family, mu, phi, MC=FALSE) {

  if (!MC) {

    if (family == "betaMode") {
      VV <- (1 + phi * mu) * (1 + phi * (1 -mu)) / (2 + phi)^2 / (3 + phi)
    }
    else if (family == "unitGompertz") {
      a <- mu^phi * (phi + 1) / phi
      E  <- a^(1/phi) * exp(a) * IGamma(1-1/phi, a)
      E2 <- a^(2/phi) * exp(a) * IGamma(1-2/phi, a)
      VV <- E2 - E^2
    }
    else if (family == "unitGamma") {
      b  <- (1+log(mu)-phi)/log(mu)
      VV <- (b / (b + 2))^phi - (b / (b + 1))^(2*phi)
    }
    else if (family == "Kumaraswamy") {
      b  <- (1 + (phi - 1) * mu^(-phi)) / phi
      E  <- b * beta(1 + 1/phi, b)
      E2 <- b * beta(1 + 2/phi, b)
      VV <- E2 - E^2
    }
  } else {

    ysim <- switch (family,
                    "betaMode" = stats::rbeta(n = 1e6, shape1 = mu * phi + 1, shape2 = (1-mu) * phi + 1),
                    "unitGompertz" = rugompertz(n = 1e6, mu = mu, phi = phi),
                    "unitGamma" = rugamma(n = 1e6, mu = mu, phi = phi),
                    "Kumaraswamy" = rkum(n = 1e6, mu = mu, phi = phi)
    )

    VV <- stats::var(ysim)
  }
  VV
}
