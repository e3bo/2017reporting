#!/usr/bin/env Rscript

library(dplyr)
source("plotting.R")

individual_sims <- readRDS("individual-sims.rds")

individual_sims %>% filter(sim %in% "2") %>% group_by() %>%
    mutate(scenario = factor(change,  levels = c("prep", "ptrans"),
                             labels = c("Reporting\nincrease",
                                        "Transmission\nincrease"))) %>%
    select(scenario, time, "Reports" = reports, "Mean" = rM,
           "Second\nfactorial\nmoment" = rSecFacMom, x) %>%
    mutate("Percent\nof change" = 100 * x) %>% select(-x) %>%
    reshape2::melt(id = c("time", "scenario")) %>%
    ggplot(aes(x = time, y = value, group = interaction(scenario, variable))) +
  geom_step() + facet_grid(variable ~ scenario, scales = "free_y") +
    scale_x_continuous(breaks = c(0, 200, 400)) +
    labs(x = "Time (weeks)", y = "Value") -> example_sims_plot

ggsave(file = "example-sims.pdf", plot = example_sims_plot, width = 84,
       height = 120, units = "mm")

sessionInfo()
