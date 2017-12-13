

n2 <- function(x){
  lambda <- x[1]
  eta <- x[2]
  T <- x[3]
  nu <- x[4]
  gamma <- (eta - lambda) * T / 2
  lambda / (nu * gamma) * (1 - (1 - exp(-2 * gamma)) / (2 * gamma))
  #(1 - (1 - exp(-2 * gamma)) / (2 * gamma))
}


library(numDeriv)
grad(n2, c(2.5, 3, 4, 1))[1]

analytic_grad <- function(x){
  lambda <- x[1]
  eta <- x[2]
  T <- x[3]
  nu <- x[4]
  gamma <- (eta - lambda) * T / 2

  df2 <- -exp((lambda - eta) * T) / (lambda - eta) - (1 - exp( (lambda - eta) * T)) * T / ( (lambda - eta)^2 * T^2)
  df1 <- eta / (eta - lambda)^2
  f1 <- lambda / (eta - lambda)
  f2 <-  (1 - (1 - exp(-2 * gamma)) / (2 * gamma))
  (f1 * df2 + df1 * f2) * 2 / T
}

analytic_grad <- function(x){
  lambda <- x[1]
  eta <- x[2]
  T <- x[3]
  nu <- x[4]
  gamma <- (eta - lambda) * T / 2

  df2 <- -exp((lambda - eta) * T) / (lambda - eta) - (1 - exp( (lambda - eta) * T)) * T / ( (lambda - eta)^2 * T^2)
  df1 <- eta / (eta - lambda)^2
  f1 <- lambda / (eta - lambda)
  f2 <-  (1 + (1 - exp( (lambda - eta) * T)) / ((lambda - eta) * T))
  (f1 * df2 + df1 * f2) * 2 / T

  A <- lambda * exp((lambda - eta) * T) / (lambda - eta)^2 + lambda * (1 - exp( (lambda - eta) * T)) * T / ( (lambda - eta)^3 * T^2)
  B <- (1 / (lambda - eta)^2 + (1 - exp( (lambda - eta) * T)) / ((lambda - eta)^3 * T)) * eta
  C <- (lambda * exp( (lambda - eta) * T) + eta) / (lambda - eta)^2 + (lambda + eta) * (1 - exp( (lambda - eta) * T)) / ((lambda - eta)^3 * T)
  C * 2 / T
}

analytic_grad( c(2.3, 2.4, .004, 1))
grad(n2, c(2.3, 2.4, .004, 1))[1]





foo <- function(k = .1, x){
(1 - exp(-k *x))/ x
}
