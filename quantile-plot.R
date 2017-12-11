#!/usr/bin/env Rscript

library(dplyr)
source("plotting.R")

quants <- readRDS("quantile-data.rds")

quants$ensemble <- factor(quants$ensemble,
                          levels = c("Individual", "Homogeneous\nensemble",
                          "Heterogeneous\nensemble"))

quants %>% mutate(scenario = factor(change,  levels = c("prep", "ptrans"),
          labels = c("Reporting\nincrease", "Transmission\nincrease"))) %>%
  mutate(stat = factor(var, levels = c("rM", "rSecFacMom"),
                       labels = c("Mean", "Second\nfactorial moment"))) %>%
    ggplot(aes(x = time, ymin = q05, ymax = q95, fill = ensemble)) +
  geom_ribbon(alpha = 0.5) +
  facet_grid(stat ~ scenario, scales = "free_y") +
  scale_x_continuous(breaks = c(200, 400)) +
  scale_fill_manual(name = "Estimate", values = pal, guide = guide_legend(title.position = "top")) +
    theme(legend.position = 'top') +
  labs(x = "Time (weeks)", y = "Value") -> quants_plot

ggsave(file = "quantiles.pdf", plot = quants_plot, width = 84, height = 100, units = "mm")

sessionInfo()
