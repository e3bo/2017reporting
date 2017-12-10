#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)

source("moment-equations.R")
source("plotting.R")

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

g <- ggplot(pdata,
            aes(x = xi, color = lambda_factor, size = lambda_factor,
                y = value))
g <- g + geom_line(lineend = "round")
g <- g + scale_color_manual(values = pal)
g <- g + scale_size_manual(values = c(2, 1.5, 1))
g <- g + facet_wrap(~moment, ncol = 1, scales = "free_y")
g <- g + theme(legend.position = "top")
leglab <- expression(paste("Trans. rate, ", lambda ~ (wk^{-1})))
g <- g + labs(x = expression(paste("Reporting probability, ", xi)),
              y = "Value", color = leglab, size = leglab)

ggsave(file = "mean-and-sfm.pdf", plot = g, height = 84, width = 84,
       units = "mm")
