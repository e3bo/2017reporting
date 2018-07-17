mean_eq <- expression(eta * T * nu * xi/ (eta - lambda))
gamma_eq <- expression((eta - lambda) * T / 2)
sfm_cases_eq <- expression(1 +
    lambda / (nu * gamma) * (1 - (1 - exp(-2 * gamma)) / (2 * gamma)))
var_nb_eq <- expression(math_mean * (1 - math_mean +
                                       (1 + 1 / phi) * (math_mean * math_sfm_cases + xi)))
cv_nb_eq <- expression(sqrt(1 / math_mean -1 + (1 + 1 / phi) * (math_sfm_cases + xi / math_mean)))
var_binom_eq <- expression(math_mean ^ 2 * math_sfm_cases + math_mean -
                             math_mean ^ 2)
cv_binom_eq <- expression(sqrt(math_sfm_cases - 1 + 1 / math_mean))
sfm_reports_nb_eq <- expression((1 + 1 / phi) *
    (math_sfm_cases + xi / math_mean))
bilinear_mom_eq <- expression(1 +
    lambda / (gamma^2 * nu) * sinh(gamma) ^ 2 * exp(-(eta - lambda) * tau))
autocor_eq <- expression((math_bmf - 1) * math_mean ^ 2 / math_var)
autocov_pref_eq <- expression(math_mean ^ 2 * lambda / (gamma^2 * nu) * sinh(gamma) ^ 2 )
eig_eq <- expression(lambda - eta)

eval_eq <- function(df, eq){
  apply(df, 1, function(x) with(as.list(x), eval(eq)))
}
