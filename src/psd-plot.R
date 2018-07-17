#!/usr/bin/env Rscript

source("moment-equations.R")
source("plotting.R")
psd_pars <- readRDS("psd_pars.rds")

sel <- !grepl("iv|lev", colnames(psd_pars)) ## avoid character types for apply()
psd_pars$math_mean <- eval_eq(psd_pars[, sel], mean_eq)
psd_pars$gamma <- eval_eq(psd_pars[, sel], gamma_eq)
psd_pars$math_sfm_cases <- eval_eq(psd_pars[, sel], sfm_cases_eq)
psd_pars$math_nb_var <- eval_eq(psd_pars[, sel], var_nb_eq)
psd_pars$math_binom_var <- eval_eq(psd_pars[, sel], var_binom_eq)
psd_pars$math_var <- ifelse(psd_pars$obsmodel > 0, psd_pars$math_binom_var,
                         psd_pars$math_nb_var)
psd_pars$math_sfm_nb_reports <- eval_eq(psd_pars[, sel], sfm_reports_nb_eq)
psd_pars$math_sfm_reports <- ifelse(psd_pars$obsmodel > 0, psd_pars$math_sfm_cases,
                                 psd_pars$math_sfm_nb_reports)
psd_pars$math_eig <- eval_eq(psd_pars[, sel], eig_eq)
psd_pars$math_ap <- eval_eq(psd_pars[, sel], autocov_pref_eq)



