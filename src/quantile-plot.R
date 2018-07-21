#!/usr/bin/env Rscript

library(dplyr)
source("plotting.R")

quants <- readRDS("quantile-data.rds")

quants$ensemble <- factor(quants$ensemble,
                          levels = c("Individual", "Homogeneous\nensemble",
                                     "Heterogeneous\nensemble"))
qplot <- function(scenario) {
  if (scenario == "single_change") {
      levs <- c("prep", "ptrans")
      labs <- c("Reporting\nincrease", "Transmission\nincrease")
  } else {
      levs <- c("prep_up_trans_up", "prep_down_trans_up")
      labs <- c("Reporting and\ntransmission\nincrease",
                "Reporting decrease,\ntransmission\nincrease")
  }
  test <- quants$change %in% levs
  quants[test, ] %>% mutate(scenario = factor(change,  levels = levs,
            labels = labs)) %>%
    mutate(stat = factor(var, levels = c("rM", "rSecFacMom"),
                         labels = c("Mean", "Second\nfactorial moment"))) %>%
      ggplot(aes(x = time, ymin = q05, ymax = q95, fill = ensemble)) +
    geom_ribbon(alpha = 0.5) +
    facet_grid(stat ~ scenario, scales = "free_y") +
    scale_x_continuous(breaks = c(200, 400)) +
      scale_fill_manual(name = "Estimate", values = pal,
                        guide = guide_legend(title.position = "top")) +
      theme(legend.position = 'top') +
      labs(x = "Time (weeks)", y = "Value") -> quants_plot
  quants_plot
}

ggsave(file = "quantiles.pdf", plot = qplot("single_change"), width = 84, height = 100,
       units = "mm")
ggsave(file = "quantiles-both-changing.pdf", plot = qplot("both"), width = 84,
       height = 100, units = "mm")

test <- quants$change %in% c("prep", "ptrans")

sessionInfo()
