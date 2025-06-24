# Functions for estimation in simulation studies ===============================

#' Estimate the ATE using the saturated, SEM, BMA-SEM, and SL-SEM estimators
#'
#' @inheritParams joint_sem
#' @param primary Name of the primary endpoint
#' @param n_boot Number of bootstrap replications to perform
#' @param ci_level The confidence level required
#' @return A data.frame with point estimates, standard error estimates, and
#'   bootstrapped percentile confidence interval bounds
#' @import sl3
#' @examples
#' \donttest{
#' data(joint_example)
#' # Using n_boot = 5 for example purposes only, use more bootstrap iterations
#' # in practice!
#' # Continuous endpoint (ATE)
#' joint_ma_sim(
#'   data0 = joint_example,
#'   endpoints = c("Y1", "Y2", "Y3_cat"),
#'   treatment = "A",
#'   primary = "Y1",
#'   n_boot = 1)
#'
#' # Binary endpoint (risk difference ATE)
#' joint_ma_sim(
#'   data0 = joint_example,
#'   endpoints = c("Y1", "Y2", "Y3_cat"),
#'   treatment = "A",
#'   primary = "Y3_cat",
#'   n_boot = 1)
#' }
#' @export
joint_ma_sim <- function(
    data0, endpoints, treatment = "A", primary, n_boot = 100,
    ci_level = 0.95, sandwich = FALSE, ...) {


  # Estimates on data0
  estimate <- estimate_ma(data0, endpoints, treatment, primary, ...)

  # Bootstrapping
  n <- nrow(data0)
  boot_results <- as.data.frame(t(replicate(
    n_boot,
    expr = {
      data_boot <- data0[sample(1:n, size = n, replace = TRUE), ]

      tryCatch(
        expr = {
          estimate_ma(data_boot, endpoints, treatment, primary, sandwich = FALSE, ...)
        },
        error = function(cond) {
          message("Error in bootstrap sampling:")
          message(conditionMessage(cond))
          out_boot_error <- rep(NA, length(estimate))
          names(out_boot_error) <- names(estimate)
          return(out_boot_error)
        }
      )
    }
  )))

  # Results
  re <- data.frame(
    Method = c("SEM", "Saturated", "SL", "BIC"),

    Estimate = estimate[c("est_sem", "est_saturated", "est_sl", "est_bic")],

    se = sqrt(c(estimate[c("v_sem", "v_saturated")], NA, NA)),

    se_boot = c(sd(boot_results$est_sem, na.rm = TRUE),
                sd(boot_results$est_saturated, na.rm = TRUE),
                sd(boot_results$est_sl, na.rm = TRUE),
                sd(boot_results$est_bic, na.rm = TRUE)),

    lb_boot = c(
      quantile(boot_results$est_sem, (1 - ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_saturated, (1 - ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_sl, (1 - ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_bic, (1 - ci_level) / 2, na.rm = TRUE)
    ),

    ub_boot = c(
      quantile(boot_results$est_sem, (1 + ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_saturated, (1 + ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_sl, (1 + ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_bic, (1 + ci_level) / 2, na.rm = TRUE)
    ),

    Omega = c(NA, NA, estimate[c("omega_sl", "omega_bic")]),

    iterations = c(
      estimate[c("nlm_iterations_sem", "nlm_iterations_saturated")],
      NA,
      NA),

    max_iterations_boot = c(
      max(c(boot_results$nlm_iterations_sem, -Inf),  na.rm = TRUE),
      max(c(boot_results$nlm_iterations_saturated, -Inf), na.rm = TRUE),
      NA,
      NA
    ),

    nlm_code = c(
      estimate[c("nlm_code_sem", "nlm_code_saturated")],
      NA,
      NA
    ),
    n_boot = c(NA, NA, sum(!is.na(boot_results$est_sl)), sum(!is.na(boot_results$est_bic)))
  )
  rownames(re) <- NULL

  return(re)
}
