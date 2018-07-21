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
  if (pars$betar > 0 && pars$betap == 0) {
    out$change <- "ptrans"
  } else if (pars$betar == 0 && pars$betap > 0) {
    out$change <- "prep"
  } else if (pars$betar > 0 && pars$betap > 0){
    out$change <- "prep_up_trans_up"
  } else if (pars$betar > 0 && pars$betap < 0){
    out$change <- "prep_down_trans_up"
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
nv_sim <- do_sims(pars, bdi, nsim = 1e4)
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

bdiuu <- bdi
pomp::coef(bdiuu)["lambda"] <- 0.5
pomp::coef(bdiuu)["betap"] <- 0.4
pomp::coef(bdiuu)["betar"] <- 0.4
pomp::coef(bdiuu)["xi"] <- 0.1
simuu <- pomp::simulate(bdiuu, as.data.frame = TRUE, times = 1:520, nsim = 1000)
simuu$change <- "prep_up_trans_up"

bdidu <- bdi
pomp::coef(bdidu)["lambda"] <- 0.5
pomp::coef(bdidu)["betap"] <- -0.4
pomp::coef(bdidu)["betar"] <- 0.4
pomp::coef(bdidu)["xi"] <- 0.9
simdu <- pomp::simulate(bdidu, as.data.frame = TRUE, times = 1:520, nsim = 1000)
simdu$change <- "prep_down_trans_up"

sim_same <- rbind(sim, simt, simuu, simdu)
saveRDS(sim_same, "sim_same.rds")

### Power spectra simulations

psd_pars <- list()
options("stringsAsFactors" = FALSE)
psd_times <- 1:520

bdi_lowrep <- bdi
pomp::coef(bdi_lowrep)["betap"] <- 0
sim_lowrep <- pomp::simulate(bdi_lowrep, as.data.frame = TRUE,
                             times = psd_times, nsim = 1000)
saveRDS(sim_lowrep, "sim_psd_lowrep.rds")
psd_pars[[1]] <- do.call(data.frame, c(list(iv = "rep", lev = "low"),
                                       as.list(pomp::coef(bdi_lowrep))))

bdi_highrep <- bdi
pomp::coef(bdi_highrep)["xi"] <- pomp::coef(bdi_highrep)["xi"] + pomp::coef(bdi)["betap"]
pomp::coef(bdi_highrep)["betap"] <- 0
sim_highrep <- pomp::simulate(bdi_highrep, as.data.frame = TRUE,
                              times = psd_times, nsim = 1000)
saveRDS(sim_highrep, "sim_psd_highrep.rds")
psd_pars[[2]] <- do.call(data.frame, c(list(iv = "rep", lev = "high"),
                                       as.list(pomp::coef(bdi_highrep))))

bdit_low <- bdit
pomp::coef(bdit_low)["betar"] <- 0
simt_low <- pomp::simulate(bdit_low, as.data.frame = TRUE,
                           times = psd_times, nsim = 1000)
saveRDS(simt_low, "sim_psd_lowtrans.rds")
psd_pars[[3]] <- do.call(data.frame, c(list(iv = "trans", lev = "low"),
                                       as.list(pomp::coef(bdit_low))))

bdit_high <- bdit
pomp::coef(bdit_high)["lambda"] <- pomp::coef(bdit_high)["lambda"] + pomp::coef(bdit)["betar"] * pomp::coef(bdit)["eta"]
pomp::coef(bdit_high)["betar"] <- 0
simt_high <- pomp::simulate(bdit_high, as.data.frame = TRUE,
                            times = psd_times, nsim = 1000)
saveRDS(simt_high, "sim_psd_hightrans.rds")
psd_pars[[4]] <- do.call(data.frame, c(list(iv = "trans", lev = "high"),
                                       as.list(pomp::coef(bdit_high))))

psd_pars <- do.call(rbind, psd_pars)
stopifnot(isTRUE(all.equal(max(abs(diff(psd_times) - psd_times[1])), 0))) ## assertion for next line
psd_pars$T <- psd_times[1]
saveRDS(psd_pars, "psd_pars.rds")

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

params3 <- params1
params3$repnum <- runif(nunits, 0.45, 0.55)
params3$xi <- runif(nunits, 0.45, 0.55)
params3$betap <- 0.4
params3$betar <- 0.4

params4 <- params1
params4$repnum <- runif(nunits, 0.45, 0.55)
params4$xi <- runif(nunits, 0.85, 0.95)
params4$betap <- -0.4
params4$betar <- 0.4

params <- rbind(params1, params2, params3, params4)
params$lambda <- params$eta * params$repnum

sim_vary <- do_sims(params, bdi)
saveRDS(sim_vary, "sim_vary.rds")

parallel::stopCluster(clust)
sessionInfo()
