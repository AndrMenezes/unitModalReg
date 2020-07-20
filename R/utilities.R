# Format output ------------------------------------------------------------

format.perc <- function(probs, digits) {
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
}

FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}



# Icomplete Upper Gamma function ------------------------------------------

IGamma <- function(a, x) stats::pgamma(x, a, lower=FALSE)*gamma(a)


# ERF function ------------------------------------------------------------

erf <- function(x) {
  2 * stats::pnorm(x * sqrt(2)) - 1
}


# Function to estimate the nuance parameter --------------------------------

est_phi <- function(family, y, yhat, p) {
  n <- length(y)
  if (family %in% c("betaMode", "GPB")) {
    sig2 <- sum((y-yhat)^2/ (n - p))
    phi  <- log(2*max(0.5*(0.25/sig2-1), 1.01)-2)
  }
  else if (family == "unitGompertz") {

    g_phi <- function(phi, y) {
      n <- length(y)
      a <- n / (sum(y^(-phi)) - n)
      n / phi + a * sum(log(y) * y ^(-phi)) - sum(log(y))
    }
    out <- tryCatch(stats::uniroot(g_phi, interval = c(1e-04, 100), y=y)[["root"]], error = function(e) NULL)
    phi <- ifelse(is.null(out), 1.5, out)
  }
  # else if (family == "unitGamma") {
  #
  #   g_phi <- function(phi, y) {
  #     n <- length(y)
  #     b <- - phi * sum(log(y)) / n
  #     n * log(b) - n * digamma(phi) + sum(log(-log(y)))
  #   }
  #   out <- tryCatch(stats::uniroot(g_phi, interval = c(1e-04, 100), y=y)[["root"]], error = function(e) NULL)
  #   phi <- ifelse(is.null(out), 1.5, out)
  # }
  else {
    phi <- 1.5
  }
  phi
}

# y <- rugompertz(1000, mu = 0.5, phi = 0.4)
# est_phi(family="unitGompertz", y=y)


# y <- rugamma(1000, mu = 0.5, phi = 1.6)
# est_phi(family="unitGamma", y=y)
