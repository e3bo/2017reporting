#!/usr/bin/env Rscript

library(ggplot2)

parsp <- readRDS("parsp.rds")

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
m1 <- reshape2::melt(parsp, id = c("parset"), measure = simvars,
                     value.name = "simulation")
m1$statistic <- gsub("^sim_", "", m1$variable)
m1$variable <- NULL
mathvars <- gsub("^sim", "math", simvars)
mathvars <- gsub("sfm$", "sfm_reports", mathvars)
m2 <- reshape2::melt(parsp, id = c("parset"), measure = mathvars,
                     value.name = "theory")
m2$statistic <- gsub("^math_", "", m2$variable)
m2$statistic <- gsub("_reports$", "", m2$statistic)
m2$variable <- NULL
m12 <- merge(m1, m2)
m12$facet <- factor(m12$statistic,
                    levels = c("mean", "sfm", "var", "bmf", "autocor"),
                    labels = c("Mean", "Second factorial moment", "Variance",
                               "Bilinear moment function", "Autocorrelation"))

g <- ggplot(m12, aes(x = theory, y = simulation)) +
      geom_abline(slope = 1, intercept = 0, col = "grey") +
        facet_wrap(~facet, ncol = 1, scales = "free") + geom_point(alpha = 0.8) +
          labs(x = "Theory", y = "Simulation")

ggsave(file = "numerical-verification.pdf", plot = g, width = 84, units = "mm")
