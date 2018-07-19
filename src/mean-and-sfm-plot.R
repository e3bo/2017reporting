#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)

source("moment-equations.R")
source("plotting.R")

pgrid <- expand.grid(xi = seq(0, 1, by = 0.001),
                     lambda = c(0.1, 0.5, 0.9))
pgrid$nu <- 1
pgrid$T <- 1
pgrid$eta <- 1
pgrid$math_mean <- eval_eq(pgrid, mean_eq)
pgrid$gamma <- eval_eq(pgrid, gamma_eq)
pgrid$math_sfm_cases <- eval_eq(pgrid, sfm_cases_eq)
pgrid$math_var <- eval_eq(pgrid, var_binom_eq)
pgrid$math_cv <- eval_eq(pgrid, cv_binom_eq)
pgrid$tau <- pgrid$T
pgrid$math_bmf <- eval_eq(pgrid, bilinear_mom_eq)
pgrid$math_autocor <- eval_eq(pgrid, autocor_eq)

pgrid %>%
  gather(math_mean, math_sfm_cases, math_var, math_cv, math_autocor,
         key = "moment", value = "value") -> pdata

pdata$moment <- factor(pdata$moment,
                       labels = c("Mean", "Second factorial moment",
                           "Variance", "Coefficient of variation",
                           "Lag-1 autocorrelation"),
                       levels = c("math_mean", "math_sfm_cases", "math_var",
                           "math_cv", "math_autocor"))
pdata$lambda_factor <- as.factor(pdata$lambda)

g <- ggplot(pdata,
            aes(x = xi, color = lambda_factor, size = lambda_factor,
                y = value))
g <- g + geom_line(lineend = "round")
g <- g + scale_color_manual(values = pal)
g <- g + scale_size_manual(values = c(2, 1.5, 1))
g <- g + facet_wrap(~moment, ncol = 1, scales = "free_y")
g <- g + theme(legend.position = "top")
leglab <- expression(paste("Trans. rate, ", lambda ~ (week^{-1})))
g <- g + labs(x = expression(paste("Reporting probability, ", xi)),
              y = "Value", color = leglab, size = leglab)

ggsave(file = "mean-and-sfm.pdf", plot = g, width = 84,
       units = "mm")

sessionInfo()
