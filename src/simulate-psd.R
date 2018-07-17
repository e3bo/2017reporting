#!/usr/bin/env Rscript

library(doParallel)

seed <- 1
set.seed(seed)

stationary_prob <- function(n, params) {
  r <- params["nu"] / params["lambda"] - 1
  ret <- n * log(params["lambda"] / params["eta"])
  ret <- ret + ((params["nu"] / params["lambda"]) *
                log((params["eta"] - params["lambda"]) / params["eta"]))
  ret <- ret + lgamma(n + r + 1)  -  lgamma(n+1) - lgamma(r + 1)
  exp(ret)
}

rstationary <- function(params) {
    rand <- runif(1)
    cum <- 0
    N <- -1
    while (cum <= rand){
        N <- N + 1
        cum <- cum + stationary_prob(N, params)
    }
    N
}

create_bdi <- function(times = 1, t0 = 0, eta = 1, lambda = 0.5, nu = 0.1,
                       betap = 0, betar = 0, N_0 = -1, phi = 1000, xi = 1,
                       obsmodel = 1,
                       covar = data.frame(x = c(0, 1),
                                          time=c(0, 520))) {
  data <- data.frame(time = times, reports = NA)
  v <- cbind(birth = c(N = 1, dct = 0, bct = 1),
             death = c(N = -1, dct = 1, bct =  0))
  rate_snippet <- pomp::Csnippet("
    if (j == 1) {
      rate = (lambda + x * betar * eta) * N + nu;
    } else {
      rate = eta * N;
    }
  ")
  rprocess <- pomp::gillespie.sim(rate.fun = rate_snippet, v = v)
  initializer <- function(params, t0, ...) {
      if (params["N_0"] < 0) {
          c(N = rstationary(params), dct = 0, bct = 0)
      } else {
          c(N = unname(params["N_0"]), dct = 0, bct = 0)
      }
  }
  params <- c(eta = eta, lambda = lambda, nu = nu, N_0 = N_0, xi = xi,
              betar = betar, betap = betap, phi = phi, obsmodel = obsmodel)
  rmeasure <- pomp::Csnippet("
    double prob = xi + x * betap;
    if (obsmodel > 0) {
      reports = rbinom(dct, prob);
    } else {
      double mean = dct * prob;
      reports = rnbinom_mu(phi, mean);
    }
  ")
  pomp::pomp(data = data, times = "time", t0 = t0,
             params = params,
             rprocess = rprocess, statenames = c("N", "dct", "bct"),
             zeronames = c("dct", "bct"),
             paramnames = names(params),
             initializer = initializer, covar = covar, tcovar = "time",
             rmeasure = rmeasure)
}

ranges <- data.frame(nu = c(0.01, 1), T = c(1, 10), xi = c(0, 1),
                     log10phi = c(-1, 2), obsmodel = c(-1, 1))
nobs <- 20
usamp <- function(x) runif(n = nobs, min = x[1], max = x[2])
pars <- data.frame(apply(ranges, 2, usamp))
pars$lambda <- seq(0.05, 0.95, by = 0.1)
pars$eta <- 1
pars$betar <- 0
pars$betap <- 0
pars$N_0 <- -1
pars$phi <- 10 ^ pars$log10phi
pars$parset <- 1:nrow(pars)

sim_with_params <- function(mod, ..., times = seq(1, 52 * 10), nsim = 1){
  pars <- list(...)
  if ("T" %in% names(pars)) times <- pars$T * c(1, 2)
  pnames <- names(pomp::coef(mod))
  pomp::coef(mod) <- unlist(pars[pnames])
  out <- pomp::simulate(mod, as.data.frame = TRUE, times = times, nsim = nsim)
  out$parset <- pars$parset
  if (pars$betar > 0 && !pars$betap > 0) {
    out$change <- "ptrans"
  } else if (!pars$betar > 0 && pars$betap > 0) {
    out$change <- "prep"
  } else {
    out$change <- "other"
  }
  out
}

clust <- parallel::makeForkCluster(nnodes = 2)
doParallel::registerDoParallel(cl = clust)
parallel::clusterSetRNGStream(cl = clust, iseed = seed)

do_sims <- function(par, md, nsim = 1){
    foreach(i = seq(1, nrow(par)), .combine = rbind) %dopar%
        do.call(sim_with_params, c(list(mod = md, nsim = nsim),
                                   as.list(par[i, ])))
}

## Simulations for numerical verification

bdi <- create_bdi(lambda = 0.5, eta = 1, N_0 = -1, nu = 1)
system.time(nv_sim <- do_sims(pars, bdi, nsim = 1e4))
saveRDS(pars, "pars.rds")
saveRDS(nv_sim, "nv_sim.rds")




## Homogeneous ensemble simulations

pomp::coef(bdi)["lambda"] <- 0.9
pomp::coef(bdi)["xi"] <- 0.1
pomp::coef(bdi)["betap"] <- 0.4
sim <- pomp::simulate(bdi, as.data.frame = TRUE, times = 1:520, nsim = 1000)
sim$change <- "prep"

bdit <- bdi
pomp::coef(bdit)["lambda"] <- 0.5
pomp::coef(bdit)["betap"] <- 0
pomp::coef(bdit)["betar"] <- 0.4
pomp::coef(bdit)["xi"] <- 0.5
simt <- pomp::simulate(bdit, as.data.frame = TRUE, times = 1:520, nsim = 1000)
simt$change <- "ptrans"

sim_same <- rbind(sim, simt)
saveRDS(sim_same, "sim_same.rds")

### Heterogeneous ensemble simulation

nunits <- 1000
params1 <- data.frame(repnum = runif(nunits, 0.85, 0.95),
                      eta = runif(nunits, 0.5 , 1.5),
                      nu = runif(nunits, 0.5, 1.5),
                      xi = runif(nunits, 0.05, 0.15),
                      betar = 0, betap = 0.4, N_0 = -1, phi = 100,
                      obsmodel = 1)
params1$parset <- 1:nrow(params1)

params2 <- params1
params2$repnum <- runif(nunits, 0.45, 0.55)
params2$xi <- runif(nunits, 0.45, 0.55)
params2$betap <- 0
params2$betar <- 0.4
params <- rbind(params1, params2)
params$lambda <- params$eta * params$repnum

sim_vary <- do_sims(params, bdi)
saveRDS(sim_vary, "sim_vary.rds")

### PSD simulations

rawspec <- function(x) {
    spec.pgram(x, taper = 0, pad = 0, fast = FALSE, demean = TRUE, detrend = FALSE, plot = FALSE)$spec
}

psd <- function(f, eig, ct, ap) {
  ct + ap * sinh(eig) / (cos(2 * pi * f) - cosh(eig))
}




bdips <- bdi
pomp::coef(bdips)["betap"] <- 0
simps <- pomp::simulate(bdips, as.data.frame = TRUE, times = 1:520, nsim = 1000)

eig <- with(as.list(pomp::coef(bdips)), lambda - eta)
mean_eq <- expression(eta * T * nu * xi/ (eta - lambda))
gamma_eq <- expression((eta - lambda) * T / 2)
sfm_cases_eq <- expression(1 + lambda / (nu * gamma) * (1 - (1 - exp(-2 * gamma)) / (2 * gamma)))
var_binom_eq <- expression(math_mean ^ 2 * math_sfm_cases + math_mean - math_mean ^ 2)
autocov_pref_eq <- expression(math_mean ^ 2 * lambda / (gamma^2 * nu) * sinh(gamma) ^ 2 )

evalf <- function(eq) with(as.list(coef(bdips)), eval(eq))
math_mean <- evalf(mean_eq)
gamma <- evalf(gamma_eq)
math_sfm_cases <- evalf(sfm_cases_eq)
var_binom <- evalf(var_binom_eq)
math_ap <- evalf(autocov_pref_eq)
math_ct <- var_binom - math_ap


anans <- psd(f = fz, eig = eig, ct = math_ct, ap = math_ap)

splt <- split(simps, simps$sim)
specs <- sapply(splt, function(x) rawspec(x$reports))
spec <- rowMeans(specs)

N <- nrow(splt[[1]])
fz <- seq(1, floor(N / 2)) / N

plot(fz, spec, log='y', type = 'l', ylab = 'Power', xlab = "Frequency (1/w)")
lines(fz, anans, col=2)

parallel::stopCluster(clust)
sessionInfo()
