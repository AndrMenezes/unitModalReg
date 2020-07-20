


# Log-likelihood of unit-Maxwell-Boltzmann --------------------------------

llumaxwell <- function(par, linkinv, X, y, sum=TRUE) {

  # location parameter (mu)
  p           <- ncol(X)
  beta        <- par[seq.int(length.out = p)]
  mu          <- linkinv(X %*% beta)
  # Log-likelihood
  theta <- (2 - mu) / log(mu)^2
  cte   <- log(1 / y) + log(2 / pi) + log(log(y)^2)
  if (sum) {
    ll <- -suppressWarnings(
      sum(cte) + 3/2 * sum(log(theta)) - 1/2 * sum(theta * log(y)^2)
    )
  } else {
    ll <- c(cte + 3/2 * log(theta) - 1/2 * theta * log(y)^2)
  }
  ll
}


# Log-likelihood of Beta in terms of mode ---------------------------------

llbetamode <- function(par, linkinv, X, y, sum=TRUE) {

  # location parameter (mu)
  p    <- ncol(X)
  beta <- par[seq.int(length.out = p)]
  mu   <- linkinv(X %*% beta)
  phi  <- exp(par[-seq.int(length.out = p)])

  # Log-likelihood
  if (sum) {
    ll <- -suppressWarnings(
      sum(stats::dbeta(x=y, shape1=mu * phi + 1, shape2=(1-mu) * phi + 1, log = TRUE))
    )
  } else {
    ll <- stats::dbeta(x=y, shape1=mu * phi + 1, shape2=(1-mu) * phi + 1, log = TRUE)
  }
  ll
}

# Log-likelihood of GPB ---------------------------------------------------

llugompertz <- function(par, linkinv, X, y, sum=TRUE) {

  # location parameter (mu)
  p     <- ncol(X)
  beta  <- par[seq.int(length.out = p)]
  mu    <- linkinv(X %*% beta)
  phi   <- par[-seq.int(length.out = p)]
  alpha <- mu^phi * (phi + 1) / phi

  # Log-likelihood
  if (sum) {
    ll <- -suppressWarnings(
      sum(log(alpha) + log(phi) - (phi + 1) * log(y) - alpha * y^(-phi) + alpha)
    )
  } else {
    ll <- c(log(alpha) + log(phi) - (phi + 1) * log(y) - alpha * y^(-phi) + alpha)
  }
  ll
}

# Log-likelihood of unit-Gamma --------------------------------------------

llugamma <- function(par, linkinv, X, y, sum=TRUE) {

  # location parameter (mu)
  p     <- ncol(X)
  beta  <- par[seq.int(length.out = p)]
  mu    <- linkinv(X %*% beta)
  phi   <- par[-seq.int(length.out = p)]
  b     <- (1+log(mu)-phi)/log(mu)

  # Log-likelihood
  if (sum) {
    ll <- -suppressWarnings(
      sum(
        phi * log(b) - lgamma(phi) + (b - 1) * log(y) + (phi - 1)*log(-log(y))
      )
    )
  } else {
    ll <- c(phi * log(b) - lgamma(phi) + (b - 1) * log(y) + (phi - 1)*log(-log(y)))
  }
  ll
}


# log-Likelihood of Kumaraswamy -------------------------------------------

llkum <- function(par, linkinv, X, y, sum=TRUE) {

  # location parameter (mu)
  p     <- ncol(X)
  beta  <- par[seq.int(length.out = p)]
  mu    <- linkinv(X %*% beta)
  phi   <- par[-seq.int(length.out = p)]

  b  <- (1 + (phi - 1) * mu^(-phi)) / phi
  # Log-likelihood
  if (sum) {
    ll <- -suppressWarnings(
      sum(
        log(phi) + log(b) + (phi - 1) * log(y) + (b - 1)*log(1-y^phi)
      )
    )
  } else {
    ll <- c(log(phi) + log(b) + (phi - 1) * log(y) + (b - 1)*log(1-y^phi))
  }
  ll
}



# Log-likelihood of GPB ---------------------------------------------------

llgpb <- function(par, linkinv, X, y, sum=TRUE) {

  # location parameter (mu)
  p    <- ncol(X)
  beta <- par[seq.int(length.out = p)]
  mu   <- linkinv(X %*% beta)
  phi  <- exp(par[-seq.int(length.out = p)])

  # Log-likelihood
  logI <- ifelse(y <= mu, log(y) - log(mu), log(1 - y) - log(1 - mu))
  cte  <- log(2*phi + 1) + log(phi + 1) - log(3*phi + 1)

  if (sum) {
    ll <- -suppressWarnings(
      sum(cte) +  phi * sum(logI) + sum(log( 2 - exp(phi*logI)))
    )
  } else {
    ll <- c(cte + phi * logI + log( 2 - exp(phi*logI)))
  }
  ll
}


