#!/usr/bin/env Rscript

library(dplyr)
library(forcats)
library(foreach)
library(ggplot2)
library(pomp)
library(reshape2)
library(tidyr)
library(zoo)

### Simulation

seed <- 1
set.seed(seed)

theme_set(theme_bw())

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
  rate_snippet <- Csnippet("
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
  rmeasure <- Csnippet("
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

pal <- c("#e2908c", "#319045", "#252525")

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
  coef(mod) <- unlist(pars[pnames])
  out <- simulate(mod, as.data.frame = TRUE, times = times, nsim = nsim)
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

second_fact_mom <- function(x){
  mean(x * (x - 1)) / mean(x) ^2
}

## Simulations for numerical verification

bdi <- create_bdi(lambda = 0.5, eta = 1, N_0 = -1, nu = 1)
nv_sim <- do_sims(pars, bdi, nsim = 1e4)

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

parallel::stopCluster(clust)



## Data processing

nv_sim %>% group_by(parset, sim) %>% summarize(m = mean(reports),
                                             m2 = mean(reports^2),
                                             prod = prod(reports)) %>%
    summarize(sim_mean = mean(m),
              sim_var = mean(m2) - mean(m) ^ 2,
              sim_sfm = (mean(m2) - mean(m)) / mean(m) ^ 2,
              sim_bmf = mean(prod) / mean(m) ^ 2,
              sim_autocor = ((mean(prod) - mean(m)^2) /
                             (mean(m2) - mean(m)^2)) ) -> sim_ests

parsp <- merge(pars, sim_ests)

mean_eq <- expression(eta * T * nu * xi/ (eta - lambda))
eval_eq <- function(df, eq){
  apply(df, 1, function(x) with(as.list(x), eval(eq)))
}
parsp$math_mean <- eval_eq(parsp, mean_eq)
gamma_eq <- expression((eta - lambda) * T / 2)
parsp$gamma <- eval_eq(parsp, gamma_eq)
sfm_cases_eq <- expression(1 +
    lambda / (nu * gamma) * (1 - (1 - exp(-2 * gamma)) / (2 * gamma)))
parsp$math_sfm_cases <- eval_eq(parsp, sfm_cases_eq)
var_nb_eq <- expression(math_mean +
    (1 + 1 / phi) * (math_mean * (math_mean * math_sfm_cases + xi)) -
      math_mean ^ 2)
parsp$math_nb_var <- eval_eq(parsp, var_nb_eq)
var_binom_eq <- expression(math_mean ^ 2 * math_sfm_cases + math_mean -
                             math_mean ^ 2)
parsp$math_binom_var <- eval_eq(parsp, var_binom_eq)
parsp$math_var <- ifelse(parsp$obsmodel > 0, parsp$math_binom_var,
                         parsp$math_nb_var)
sfm_reports_nb_eq <- expression((1 + 1 / phi) *
    (math_mean * (math_mean * math_sfm_cases + xi)) / math_mean ^ 2)
parsp$math_sfm_nb_reports <- eval_eq(parsp, sfm_reports_nb_eq)
parsp$math_sfm_reports <- ifelse(parsp$obsmodel > 0, parsp$math_sfm_cases,
                                 parsp$math_sfm_nb_reports)
bilinear_mom_eq <- expression(1 +
    lambda / (gamma^2 * nu) * sinh(gamma) ^ 2 * exp(-(eta - lambda) * tau))
parsp$tau <- parsp$T
parsp$math_bmf <- eval_eq(parsp, bilinear_mom_eq)
autocor_eq <- expression((math_bmf - 1) * math_mean ^ 2 / math_var)
parsp$math_autocor <- eval_eq(parsp, autocor_eq)

## Plotting

simvars <- colnames(parsp)[grep("^sim", colnames(parsp))]
m1 <- melt(parsp, id = c("parset"), measure = simvars,
           value.name = "simulation")
m1$statistic <- gsub("^sim_", "", m1$variable)
m1$variable <- NULL
mathvars <- gsub("^sim", "math", simvars)
mathvars <- gsub("sfm$", "sfm_reports", mathvars)
m2 <- melt(parsp, id = c("parset"), measure = mathvars, value.name = "theory")
m2$statistic <- gsub("^math_", "", m2$variable)
m2$statistic <- gsub("_reports$", "", m2$statistic)
m2$variable <- NULL
m12 <- merge(m1, m2)
m12$facet <- factor(m12$statistic,
                    levels = c("mean", "sfm", "var", "bmf", "autocor"),
                    labels = c("Mean", "Second factorial moment", "Variance",
                               "Bilinear moment function", "Autocorrelation"))


g <- ggplot(m12, aes(x = theory, y = simulation)) +
#  scale_x_continuous(trans = "log1p") +
#    scale_y_continuous(trans = "log1p") +
      geom_abline(slope = 1, intercept = 0, col = "grey") +
        facet_wrap(~facet, ncol = 1, scales = "free") + geom_point(alpha = 0.8) +
          labs(x = "Theory", y = "Simulation")

ggsave(file = "numerical-verification.pdf", plot = g, width = 84, units = "mm")

x_eq <- expression(eta * T * nu * xi)
parsp$mean_xcoord <- eval_eq(parsp, x_eq)
mean_group_eq <- expression(eta - lambda)
parsp$mean_group <- eval_eq(parsp, mean_group_eq)

gl <- levels(as.factor(parsp$mean_group))
parsp$"reporting\nmodel" <- ifelse(parsp$obsmodel > 0,
                                   "binomial", "negative\nbinomial")
ggplot(parsp, aes(x = mean_xcoord, y = math_mean,
                  color = as.factor(mean_group))) +
    geom_point(aes(y = sim_mean, shape = `reporting\nmodel`)) +
    scale_x_log10(limits = c(1e-2, 1e1), expand = c(0, 0)) +
    scale_y_log10(limits = c(1e-2, 1e2), expand = c(0, 0)) +
    labs(x = x_eq) + labs(y = "Mean no. reports", color = mean_group_eq) -> gg

for (lev in gl){
    eval(substitute(gg <- gg +  stat_function(fun = function(x, el) x / el,
                                              args = list(el = a),
                                              aes(color = b),
                                              inherit.aes = FALSE),
                    list(a = as.numeric(lev), b = lev)))
}
ggsave(file = "mean.pdf", plot = gg)

## By using substitute, we avoid having to write out calls explicit calls
## to get the line colors to match those of the points


mean_plot_fun <- function(x, ...){
  params <- c(list(xi=x), list(...))
  with(params, eval(mean_eq))
}

pgrid <- expand.grid(xi = seq(0, 1, by = 0.1),
                     lambda = c(0.1, 0.5, 0.9))
pgrid$nu <- 1
pgrid$T <- 1
pgrid$eta <- 1
pgrid$mean <- eval_eq(pgrid, mean_eq)
pgrid$gamma <- eval_eq(pgrid, gamma_eq)
pgrid$second_factorial_mom <- eval_eq(pgrid, sfm_cases_eq)
pgrid %>%
  gather(mean, second_factorial_mom, key = "moment", value = "value") -> pdata

pdata$moment <- factor(pdata$moment,
                       labels = c("Mean", "Second factorial moment"),
                       levels = c("mean", "second_factorial_mom"))
pdata$lambda_factor <- as.factor(pdata$lambda)
#pdata$lambda_factor) <- factor(pdata$lambda, labels = rev(labels(levels(pdata$lambda_factor))

g <- ggplot(pdata, aes(x = xi, color = lambda_factor, size = lambda_factor, y = value))
g <- g + geom_line(lineend = "round")
g <- g + scale_color_manual(values = pal)
g <- g + scale_size_manual(values = c(2, 1.5, 1))
g <- g + facet_wrap(~moment, ncol = 1, scales = "free_y")
g <- g + theme(legend.position = "top")
leglab <- expression(paste("Trans. rate, ", lambda ~ (wk^{-1})))
g <- g + labs(x = expression(paste("Reporting probability, ", xi)),
              y = "Value", color = leglab, size = leglab)

ggsave(file = "mean-and-sfm.pdf", plot = g, height = 84, width = 84, units = "mm")

### Data processing

win_size <- 52
sim_same %>% group_by(change, sim) %>%
  mutate(rM = rollmean(reports, k = win_size,
             fill = NA, align = "right")) %>%
    mutate(rV = rollapplyr(reports, width = win_size,
               FUN = var, fill = NA)) %>%
      mutate(rSecFacMom = rollapplyr(reports, width = win_size,
                 FUN = second_fact_mom, fill = NA)) %>%
  mutate(scenario = factor(change,  levels = c("prep", "ptrans"),
             labels = c("Reporting\nincrease",
                 "Transmission\nincrease"))) -> simb2





simb2 %>% filter(time > win_size - 0.5) %>% group_by(change, time) %>%
  dplyr::summarize(mean_rM = mean(rM),
                   mean_rsfm = mean(rSecFacMom),
                   q05_rM = quantile(rM, prob = 0.05),
                   q95_rM = quantile(rM, prob = 0.95),
                   q05_rsfm = quantile(rSecFacMom, prob = 0.05),
                   q95_rsfm = quantile(rSecFacMom, prob = 0.95)) -> foo

foo %>% select(change, time, q05_rM, q05_rsfm) %>%
  gather(q05_rM, q05_rsfm, key = "var", value = "q05") %>%
    mutate(var = stringr::str_replace(var, "q05_", "")) -> foo1

foo %>% select(change, time, q95_rM, q95_rsfm) %>%
  gather(q95_rM, q95_rsfm, key = "var", value = "q95") %>%
    mutate(var = stringr::str_replace(var, "q95_", "")) -> foo2

foom <- merge(foo1, foo2)

foom %>% mutate(scenario = factor(change,  levels = c("prep", "ptrans"),
      labels = c("Reporting\nincrease", "Transmission\nincrease"))) %>%
  mutate(stat = factor(var, levels = c("rM", "rsfm"),
             labels = c("Mean", "Second\nfactorial moment"))) -> iquants


## Plotting

simb2 %>% filter(sim %in% c("2")) %>% group_by() %>%
  select(scenario, time, "Reports" = reports,
         "Second\nfactorial\nmoment" = rSecFacMom, "Mean" = rM,
         x) %>% mutate("Percent\nof change" = 100 * x) %>% select(-x) %>%
    melt(id = c("time", "scenario")) -> simb3




simb3 %>% ggplot(aes(x = time, y = value,
                     group = interaction(scenario, variable))) +
  geom_step() + facet_grid(variable ~ scenario, scales = "free_y") +
    scale_x_continuous(breaks = c(0, 200, 400)) +
    labs(x = "Time (weeks)", y = "Value") -> example_sims_plot

ggsave(file = "example-sims.pdf", plot = example_sims_plot, width = 84,
       height = 120, units = "mm")

## Data processing

samp_stat <- function(x, w) {
  colMeans(x[w, ], na.rm = TRUE)
}

boot_fn <- function(df, R = 300) {
  voi <- c("rM", "rSecFacMom")
  data <- df[, voi, drop = FALSE]
  bdf <- boot::boot(data = data, statistic = samp_stat, R = R)
  ret1 <- apply(bdf$t, 2, quantile, prob = 0.95)
  ret2 <- apply(bdf$t, 2, quantile, prob = 0.05)
  ret <- data.frame(q05 = ret2, q95 = ret1)
  ret$var <- voi
  ret
}

simb2 %>% filter(time > 52) %>% group_by(change, time) %>%
  do(boot_fn(.)) -> bquants

bquants %>% mutate(scenario = factor(change,  levels = c("prep", "ptrans"),
          labels = c("Reporting\nincrease", "Transmission\nincrease"))) %>%
  mutate(stat = factor(var, levels = c("rM", "rSecFacMom"),
             labels = c("Mean", "Second\nfactorial moment"))) -> bquants

bquants$ensemble <- "Homogeneous\nensemble"
iquants$ensemble <- "Individual"
quants <- merge(bquants, iquants, all = TRUE)


## Data processing

win_size <- 52
out %>% group_by(change, parset) %>%
  mutate(rM=rollmean(reports, k = win_size, fill = NA, align = "right")) %>%
    mutate(rV=rollapplyr(reports, width = win_size, FUN = var, fill = NA)) %>%
      mutate(rSecFacMom = rollapplyr(reports, width = win_size,
                 FUN = second_fact_mom, fill = NA)) -> out2

out2 %>% filter(time > 52) %>% group_by(change, time) %>%
  do(boot_fn(.)) -> outquants

outquants %>% mutate(scenario = factor(change,  levels = c("prep", "ptrans"),
  labels = c("Reporting\nincrease", "Transmission\nincrease"))) %>%
  mutate(stat = factor(var, levels = c("rM", "rSecFacMom"),
             labels = c("Mean", "Second\nfactorial moment"))) -> oquants
oquants$ensemble <- "Heterogeneous\nensemble"
qq <- merge(quants, oquants, all = TRUE)
qq$ensemble <- factor(qq$ensemble,
                      levels = c("Individual", "Homogeneous\nensemble",
                          "Heterogeneous\nensemble"))

## Plotting

ggplot(qq, aes(x = time, ymin = q05, ymax = q95, fill = ensemble)) +
  geom_ribbon(alpha = 0.5) +
  facet_grid(stat ~ scenario, scales = "free_y") +
  scale_x_continuous(breaks = c(200, 400)) +
  scale_fill_manual(name = "Estimate", values = pal, guide = guide_legend(title.position = "top")) +
    theme(legend.position = 'top') +
#      guide_legend(title.position = "top") +
  labs(x = "Time (weeks)", y = "Value") -> quants_plot

ggsave(file = "quantiles.pdf", plot = quants_plot, width = 84, height = 100, units = "mm")



## Reproducibility

save.image(file="reporting-vs-transmission.RData")

sessionInfo()

