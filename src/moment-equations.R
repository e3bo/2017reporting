mean_eq <- expression(eta * T * nu * xi/ (eta - lambda))
gamma_eq <- expression((eta - lambda) * T / 2)
sfm_cases_eq <- expression(1 +
    lambda / (nu * gamma) * (1 - (1 - exp(-2 * gamma)) / (2 * gamma)))
var_nb_eq <- expression(math_mean +
    (1 + 1 / phi) * (math_mean * (math_mean * math_sfm_cases + xi)) -
      math_mean ^ 2)
var_binom_eq <- expression(math_mean ^ 2 * math_sfm_cases + math_mean -
                             math_mean ^ 2)
sfm_reports_nb_eq <- expression((1 + 1 / phi) *
    (math_sfm_cases + xi / math_mean))
bilinear_mom_eq <- expression(1 +
    lambda / (gamma^2 * nu) * sinh(gamma) ^ 2 * exp(-(eta - lambda) * tau))
autocor_eq <- expression((math_bmf - 1) * math_mean ^ 2 / math_var)

eval_eq <- function(df, eq){
  apply(df, 1, function(x) with(as.list(x), eval(eq)))
}
