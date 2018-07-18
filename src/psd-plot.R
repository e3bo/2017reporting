#!/usr/bin/env Rscript

source("moment-equations.R")
source("plotting.R")
psd_pars <- readRDS("psd_pars.rds")

dsel <- which(grepl("iv|lev", colnames(psd_pars))) ## avoid character for apply()
psd_pars$math_mean <- eval_eq(psd_pars[, -dsel], mean_eq)
psd_pars$gamma <- eval_eq(psd_pars[, -dsel], gamma_eq)
psd_pars$math_sfm_cases <- eval_eq(psd_pars[, -dsel], sfm_cases_eq)
psd_pars$math_nb_var <- eval_eq(psd_pars[, -dsel], var_nb_eq)
psd_pars$math_binom_var <- eval_eq(psd_pars[, -dsel], var_binom_eq)
psd_pars$math_var <- ifelse(psd_pars$obsmodel > 0, psd_pars$math_binom_var,
                         psd_pars$math_nb_var)
psd_pars$math_sfm_nb_reports <- eval_eq(psd_pars[, -dsel], sfm_reports_nb_eq)
psd_pars$math_sfm_reports <- ifelse(psd_pars$obsmodel > 0,
                                    psd_pars$math_sfm_cases,
                                 psd_pars$math_sfm_nb_reports)
psd_pars$math_eig <- eval_eq(psd_pars[, -dsel], eig_eq)
psd_pars$math_ap <- eval_eq(psd_pars[, -dsel], autocov_pref_eq)

## Plotting

rawspec <- function(x) {
  spec.pgram(x, taper = 0, pad = 0, fast = FALSE, demean = TRUE,
             detrend = FALSE, plot = FALSE)$spec
}

psd <- function(f, eig, var, ap) {
  ct <- var - ap
  ct + ap * sinh(eig) / (cos(2 * pi * f) - cosh(eig))
}

spec_df <- function(pars){
  simfile <- paste0("sim_psd_", pars$lev, pars$iv, ".rds")
  sim <- readRDS(simfile)
  splt <- split(sim, sim$sim)
  specs <- sapply(splt, function(x) rawspec(x$reports))
  spec <- rowMeans(specs)
  N <- nrow(splt[[1]])
  fz <- seq(1, floor(N / 2)) / N
  aspec <- psd(f = fz, eig = pars$math_eig, var = pars$math_var,
               ap = pars$math_ap)
  data.frame(fz = fz, spec_num = spec, spec_an = aspec, iv = pars$iv,
             lev = pars$lev)
}

sdf <- list()
for(i in seq_len(nrow(psd_pars))){
  sdf[[i]] <- spec_df(pars = psd_pars[i, ])
}

sdata <- do.call(rbind, sdf)
sdm <- reshape2::melt(sdata, id.vars = c("iv", "lev", "fz"),
                      value.name = "Power")

sdm$Method <- factor(sdm$variable, levels = c("spec_num", "spec_an"),
                     labels = c("Numeric", "Analytic"))
sdm$iv_facet <- factor(sdm$iv, levels = c("rep", "trans"),
                       labels = c("Reporting\nincrease",
                           "Transmission\nincrease"))
sdm$lev_facet <- factor(sdm$lev, levels = c("low", "high"),
                        labels = c("Before", "After"))

g <- ggplot(data = sdm, aes(x = fz, y = Power, color = Method))
g <- g + geom_line() + facet_grid(lev_facet ~ iv_facet, scales = "free_y")
g <- g + scale_y_log10() + xlab("Frequency (1 / week)")
g <- g + scale_color_manual(values = pal[-2],
                            guide = guide_legend(title.position = "left"))
g <- g + theme(legend.position = 'top')

ggsave(file = "psd.pdf", plot = g, width = 84, height = 100, units = "mm")
