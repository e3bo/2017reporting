#!/usr/bin/env Rscript

library(dplyr)

nv_sim <- readRDS("nv_sim.rds")
pars <- readRDS("pars.rds")

nv_sim %>% group_by(parset, sim) %>% summarize(m = mean(reports),
                                             m2 = mean(reports^2),
                                             prod = prod(reports)) %>%
    summarize(sim_mean = mean(m),
              sim_var = mean(m2) - mean(m) ^ 2,
              sim_cv = sqrt(mean(m2) - mean(m) ^ 2) / mean(m),
              sim_sfm = (mean(m2) - mean(m)) / mean(m) ^ 2,
              sim_bmf = mean(prod) / mean(m) ^ 2,
              sim_autocor = ((mean(prod) - mean(m)^2) /
                             (mean(m2) - mean(m)^2)) ) -> sim_ests

parsp <- merge(pars, sim_ests)
saveRDS(parsp, "parsp.rds")

sessionInfo()
