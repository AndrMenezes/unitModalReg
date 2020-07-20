#' @name umaxwell
#' @aliases umaxwell dumaxwell pumaxwell qumaxwell rumaxwell
#'
#' @title unit-maxwell distribution
#'
#' @description Density function, distribution function, quantile function, random number generation function
#' for the unit-maxwell distribution re-parametrized in terms of the mode.
#'
#' @author André Felipe B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the mode.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#'
#' @return \code{dumaxwell} gives the density, \code{pumaxwell} gives the distribution function,
#' \code{qumaxwell} gives the quantile function and \code{rumaxwell} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @examples
#' set.seed(6969)
#' x <- rumaxwell(n = 1e4, mu = 0.5)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.01)
#' hist(x, prob = TRUE, main = 'unit-maxwell distribution')
#' lines(S, dumaxwell(x = S, mu = 0.5), col = 2)
#' plot(ecdf(x))
#' lines(S, pumaxwell(q = S, mu = 0.5), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qumaxwell(p = S, mu = 0.5), col = 2)



# Density -----------------------------------------------------------------
#' @rdname umaxwell
#' @export
dumaxwell <- function(x, mu, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1)
  theta <- (2 - log(mu)) / (log(mu))^2
  pdf <- 1 / x * sqrt(2 / pi) * theta^(3/2) * log(x)^2 * exp(-1/2 * theta * log(x)^2 )
  if (log) {
    log(pdf)
  } else {
    pdf
  }
}

# CDF ---------------------------------------------------------------------
#' @rdname umaxwell
#' @export
pumaxwell <- function(q, mu, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(q >= 0, q <= 1, mu > 0, mu < 1)
  theta <- (2 - log(mu)) / (log(mu))^2
  p <- 1 - (erf(-log(q) * sqrt(theta / 2) ) + log(q) * sqrt(2 * theta / pi) * exp(-1/2 * theta * log(q)^2))
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

# Quantile ----------------------------------------------------------------
#' @rdname umaxwell
#' @export
qumaxwell <- function(p, mu, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p >= 0, p <= 1, mu > 0, mu < 1)
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- exp(p)
  }
  theta <- (2 - log(mu)) / (log(mu))^2
  exp(-sqrt(2 / theta * stats::qgamma(1-p, shape=3/2, scale=1)))
}

# Random numbers ----------------------------------------------------------
#' @rdname umaxwell
#' @export
rumaxwell <- function(n, mu)
{
  stopifnot(n > 0, mu > 0, mu < 1)
  u <- stats::runif(n)
  qumaxwell(u, mu)
}


#' @name ugompertz
#' @aliases ugompertz dugompertz pugompertz qugompertz rugompertz
#'
#' @title unit-Gompertz distribution
#'
#' @description Density function, distribution function, quantile function, random number generation function
#' for the unit-Gompertz distribution re-parametrized in terms of the mode.
#'
#' @author André Felipe B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the mode.
#' @param phi shape parameter.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#'
#' @return \code{dugompertz} gives the density, \code{pugompertz} gives the distribution function,
#' \code{qugompertz} gives the quantile function and \code{rugompertz} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @examples
#' set.seed(6969)
#' x <- rugompertz(n = 1e4, mu = 0.5, phi = 1.0)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.01)
#' hist(x, prob = TRUE, main = 'unit-gompertz distribution')
#' lines(S, dugompertz(x = S, mu = 0.5, phi = 1.0), col = 2)
#' plot(ecdf(x))
#' lines(S, pugompertz(q = S, mu = 0.5, phi = 1.0), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qugompertz(p = S, mu = 0.5, phi = 1.0), col = 2)



# Density -----------------------------------------------------------------
#' @rdname ugompertz
#' @export
dugompertz <- function(x, mu, phi, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, phi > 0)
  alpha <- mu^phi * (phi + 1) / phi
  pdf   <- alpha * phi * x^(-phi - 1) * exp(- alpha * (x^(-phi) - 1))
  if (log) {
    log(pdf)
  } else {
    pdf
  }
}

# CDF ---------------------------------------------------------------------
#' @rdname ugompertz
#' @export
pugompertz <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(q >= 0, q <= 1, mu > 0, mu < 1)
  alpha <- mu^phi * (phi + 1) / phi
  p <- exp(-alpha * (q ^ (-phi) - 1))
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

# Quantile ----------------------------------------------------------------
#' @rdname ugompertz
#' @export
qugompertz <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p >= 0, p <= 1, mu > 0, mu < 1, phi > 0)
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- exp(p)
  }
  alpha <- mu^phi * (phi + 1) / phi
  exp(-(log(alpha - log(p)) - log(alpha)) / phi)
}

# Random numbers ----------------------------------------------------------
#' @rdname ugompertz
#' @export
rugompertz <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0)
  u <- stats::runif(n)
  qugompertz(u, mu, phi)
}


# unit-Gamma --------------------------------------------------------------

#' @name ugamma
#' @aliases ugamma dugamma pugamma qugamma rugamma
#'
#' @title unit-Gamma distribution
#'
#' @description Density function, distribution function, quantile function, random number generation function
#' for the unit-Gamma distribution re-parametrized in terms of the mode.
#'
#' @author André Felipe B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the mode.
#' @param phi shape parameter.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#'
#' @return \code{dugamma} gives the density, \code{pugamma} gives the distribution function,
#' \code{qugamma} gives the quantile function and \code{rugamma} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @examples
#' set.seed(6969)
#' x <- rugamma(n = 1e4, mu = 0.5, phi = 3.0)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.01)
#' hist(x, prob = TRUE, main = 'unit-gamma distribution')
#' lines(S, dugamma(x = S, mu = 0.5, phi = 3.0), col = 2)
#' plot(ecdf(x))
#' lines(S, pugamma(q = S, mu = 0.5, phi = 3.0), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qugamma(p = S, mu = 0.5, phi = 3.0), col = 2)



# Density -----------------------------------------------------------------
#' @rdname ugamma
#' @export
dugamma <- function(x, mu, phi, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, phi > 0)
  b   <- (1+log(mu)-phi)/log(mu)
  pdf <- b ^ phi * x ^ (b - 1) * (-log(x)) ^ (phi - 1) / gamma(phi)
  if (log) {
    log(pdf)
  } else {
    pdf
  }
}

# CDF ---------------------------------------------------------------------
#' @rdname ugamma
#' @export
pugamma <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(q >= 0, q <= 1, mu > 0, mu < 1)
  b <- (1+log(mu)-phi)/log(mu)
  p <- stats::pgamma(q = -log(q), shape = phi, rate = b, lower.tail = F)
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

# Quantile ----------------------------------------------------------------
#' @rdname ugamma
#' @export
qugamma <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p >= 0, p <= 1, mu > 0, mu < 1, phi > 0)
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- exp(p)
  }
  b <- (1+log(mu)-phi)/log(mu)
  exp(-stats::qgamma(p = p, shape = phi, rate = b, lower.tail = F))
}

# Random numbers ----------------------------------------------------------
#' @rdname ugamma
#' @export
rugamma <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0)
  u <- stats::runif(n)
  qugamma(u, mu, phi)
}


# Kumaraswamy -------------------------------------------------------------

#' @name Kumaraswamy
#' @aliases kum dkum pkum qkum rkum
#'
#' @title Kumaraswamy distribution
#'
#' @description Density function, distribution function, quantile function, random number generation function
#' for the Kumaraswamy distribution re-parametrized in terms of the mode.
#'
#' @author André Felipe B. Menezes \email{andrefelipemaringa@gmail.com}
#'
#'
#' @param x,q vector of positive quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu location parameter indicating the mode.
#' @param phi shape parameter.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)}.
#'
#'
#' @return \code{dkum} gives the density, \code{pkum} gives the distribution function,
#' \code{qkum} gives the quantile function and \code{rkum} generates random deviates.
#'
#' @return Invalid arguments will return an error message.
#'
#' @examples
#' set.seed(6969)
#' x <- rkum(n = 1e4, mu = 0.5, phi = 3.0)
#' R <- range(x)
#' S <- seq(from = R[1], to = R[2], by = 0.01)
#' hist(x, prob = TRUE, main = 'unit-gamma distribution')
#' lines(S, dkum(x = S, mu = 0.5, phi = 3.0), col = 2)
#' plot(ecdf(x))
#' lines(S, pkum(q = S, mu = 0.5, phi = 3.0), col = 2)
#' plot(quantile(x, probs = S), type = "l")
#' lines(qkum(p = S, mu = 0.5, phi = 3.0), col = 2)



# Density -----------------------------------------------------------------
#' @rdname Kumaraswamy
#' @export
dkum <- function(x, mu, phi, log = FALSE) {
  stopifnot(x > 0, x < 1, mu > 0, mu < 1, phi > 0)
  b  <- (1 + (phi - 1) * mu^(-phi)) / phi
  pdf <- phi * b * x ^ (phi - 1) * (1 - x ^ phi) ^ (b - 1)
  if (log) {
    log(pdf)
  } else {
    pdf
  }
}

# CDF ---------------------------------------------------------------------
#' @rdname Kumaraswamy
#' @export
pkum <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(q >= 0, q <= 1, mu > 0, mu < 1)
  b <- (1 + (phi - 1) * mu^(-phi)) / phi
  p <- 1 - (1 - q ^ phi) ^ b
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

# Quantile ----------------------------------------------------------------
#' @rdname Kumaraswamy
#' @export
qkum <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(p >= 0, p <= 1, mu > 0, mu < 1, phi > 0)
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- exp(p)
  }
  b  <- (1 + (phi - 1) * mu^(-phi)) / phi
  (1 - (1 - p) ^ (1 / b)) ^ (1 / phi)
}

# Random numbers ----------------------------------------------------------
#' @rdname Kumaraswamy
#' @export
rkum <- function(n, mu, phi) {
  stopifnot(n > 0, mu > 0, mu < 1, phi > 0)
  u <- stats::runif(n)
  qkum(u, mu, phi)
}


