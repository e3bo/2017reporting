#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(zoo)

sim_same <- readRDS("sim_same.rds")
sim_vary <- readRDS("sim_vary.rds")

second_fact_mom <- function(x){
  mean(x * (x - 1)) / mean(x) ^2
}

win_size <- 52
sim_same %>% group_by(change, sim) %>%
  mutate(rM = rollmean(reports, k = win_size,
             fill = NA, align = "right")) %>%
    mutate(rV = rollapplyr(reports, width = win_size,
               FUN = var, fill = NA)) %>%
      mutate(rSecFacMom = rollapplyr(reports, width = win_size,
                 FUN = second_fact_mom, fill = NA)) -> mwe_same

saveRDS(mwe_same, "mwe-same.rds")

mwe_same %>% filter(time > win_size - 0.5) %>%
    group_by(change, time) %>%
  dplyr::summarize(mean_rM = mean(rM, na.rm = TRUE),
                   mean_rSecFacMom = mean(rSecFacMom, na.rm = TRUE),
                   q05_rM = quantile(rM, prob = 0.05, na.rm = TRUE),
                   q95_rM = quantile(rM, prob = 0.95, na.rm = TRUE),
                   q05_rSecFacMom = quantile(rSecFacMom, prob = 0.05, na.rm = TRUE),
                   q95_rSecFacMom = quantile(rSecFacMom, prob = 0.95, na.rm = TRUE)) -> foo

foo %>% select(change, time, q05_rM, q05_rSecFacMom) %>%
  gather(q05_rM, q05_rSecFacMom, key = "var", value = "q05") %>%
    mutate(var = stringr::str_replace(var, "q05_", "")) -> foo1

foo %>% select(change, time, q95_rM, q95_rSecFacMom) %>%
  gather(q95_rM, q95_rSecFacMom, key = "var", value = "q95") %>%
    mutate(var = stringr::str_replace(var, "q95_", "")) -> foo2

foom <- merge(foo1, foo2)
foom -> iquants

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

mwe_same %>% filter(time > 52) %>% group_by(change, time) %>%
  do(boot_fn(.)) -> bquants

sim_vary %>% group_by(change, parset) %>%
  mutate(rM=rollmean(reports, k = win_size, fill = NA, align = "right")) %>%
    mutate(rV=rollapplyr(reports, width = win_size, FUN = var, fill = NA)) %>%
      mutate(rSecFacMom = rollapplyr(reports, width = win_size,
                 FUN = second_fact_mom, fill = NA)) -> mwe_vary

mwe_vary %>% filter(time > 52) %>% group_by(change, time) %>%
  do(boot_fn(.)) -> vquants

vquants$ensemble <- "Heterogeneous\nensemble"
bquants$ensemble <- "Homogeneous\nensemble"
iquants$ensemble <- "Individual"
quants <- merge(bquants, iquants, all = TRUE)
quants <- merge(quants, vquants, all = TRUE)

saveRDS(quants, "quantile-data.rds")

sessionInfo()
