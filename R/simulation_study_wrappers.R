# Functions for estimation in simulation studies
#
# Fit multiple models for one data set and report the results under each model.
# In practice, the model and estimation strategy should be defined in advance
# and we do not recommend using these functions outside of comparative
# simulation studies.

#' Estimate the ATE with various methods for simulation studies
#'
#' @inheritParams joint_sem
#' @param n_boot Number of bootstrap replications to perform
#' @param ci_level The confidence level required
#' @param categorical Character vector of endpoint names to be treated as
#'   categorical (binary or ordinal) variables
#' @return A data.frame with point estimates, standard error estimates, test
#'   statistics, and 95% confidence interval bounds
#' @importFrom stats pnorm uniroot var quantile qnorm
#' @examples
#' data(joint_example)
#' # Using n_boot = 5 for example purposes only, use more bootstrap iterations
#' # in practice!
#' \donttest{
#'   joint_ma_sim(joint_example, c("Y1", "Y2", "Y3"), n_boot = 5)
#' }
#'
#' @export
joint_ma_sim <- function(
    data0, endpoints, categorical = c(), treatment = "A", n_boot = 100,
    ci_level = 0.95, sandwich = FALSE, ...) {

  est_ma <- estimate_ma(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, sandwich = sandwich, ...)

  # Estimate sampling distribution of ma estimator treating weights as fixed
  v_ma_fixed <- drop(est_ma$omega %*% (est_ma$ate_ma^2 + est_ma$v_models)) -
    est_ma$ate_ma^2

  Z_ma_fixed <- est_ma$ate_ma / sqrt(v_ma_fixed)

  alpha <- 1 - ci_level

  # Estimated CDF
  F_ma <- function(x) {
    drop(est_ma$omega %*%
           pnorm(x, mean = est_ma$ate_models, sd = sqrt(est_ma$v_models))
    )
  }

  # Estimate %tile CI from CDF
  lb_fixed <- tryCatch(
    uniroot(
      f = function(x) F_ma(x) - alpha / 2,
      interval = est_ma$ate_ma + c(-5, 5) * sqrt(v_ma_fixed),
      extendInt = "upX"
    )$root,
    error = function(err) NA
  )

  ub_fixed <- tryCatch(
    uniroot(
      f = function(x) F_ma(x) - (1 - alpha / 2),
      interval = est_ma$ate_ma + c(-5, 5) * sqrt(v_ma_fixed),
      extendInt = "upX"
    )$root,
    error = function(err) NA
  )


  # Estimate sampling distribution via bootstrap

  boot_ests <- replicate(
    n_boot, do_boot_ma(data0, endpoints, categorical, treatment, ...)$ate_ma
  )
  v_ma_boot <- var(boot_ests)
  Z_ma_boot <- est_ma$ate_ma / sqrt(v_ma_boot)
  lb_boot <- quantile(boot_ests, alpha / 2)
  ub_boot <- quantile(boot_ests, 1 - alpha / 2)


  # Return point estimates, standard error estimates, test statistics, and
  # 1-alpha% CI bounds
  re <- data.frame(
    Method = c("SEM", "Saturated", "MA Fixed", "MA Boot", "MA Boot %Tile"),
    estimate = c(est_ma$ate_models, rep(est_ma$ate_ma, 3)),
    se = sqrt(c(est_ma$v_models, v_ma_fixed, rep(v_ma_boot, 2)))
  )
  re$Z <- re$estimate / re$se
  re$lb <- c(
    qnorm(alpha / 2, mean = est_ma$ate_models, sd = sqrt(est_ma$v_models)),
    lb_fixed,
    qnorm(alpha / 2, mean = est_ma$ate_ma, sd = sqrt(v_ma_boot)),
    lb_boot
  )
  re$ub <- c(
    qnorm(1 - alpha / 2, mean = est_ma$ate_models, sd = sqrt(est_ma$v_models)),
    ub_fixed,
    qnorm(1 - alpha / 2, mean = est_ma$ate_ma, sd = sqrt(v_ma_boot)),
    ub_boot
  )
  re$ll <- c(est_ma$lls, rep(NA, nrow(re) - 2))
  return(re)
}

#' Estimate the ATE with various methods for simulation studies
#'
#' @inheritParams joint_sem
#' @inheritParams estimate_ma_2
#' @param categorical Character vector of endpoint names to be treated as
#'   categorical (binary or ordinal) variables
#' @return A data.frame with point estimates for various methods
#' @examples
#' \donttest{
#' data(joint_example)
#' joint_ma_sim_2(
#'     joint_example,
#'     endpoints = c("Y1", "Y2", "Y3"),
#'     treatment = "A",
#'     primary = "Y1"
#'   )
#' }
#'
#' @export
joint_ma_sim_2 <- function(
    data0, endpoints, categorical = c(), treatment = "A", primary, sandwich = FALSE,
    ...) {


  estimate <- estimate_ma_2(data0, endpoints, categorical, treatment, primary, sandwich, ...)

  re <- data.frame(
    Method = c("SEM", "Saturated", "MA BIC", "MA SL"),
    estimate = c(estimate$ate_models, estimate$ate_ma_bic, estimate$ate_ma_sl),
    omega_bic = c(estimate$omega_bic, NA, NA),
    omega_sl  = c(estimate$omega_sl, NA, NA),
    ll = c(estimate$lls, NA, NA)
  )

  return(re)
}

#' Estimate the ATE with various methods for simulation studies
#'
#' @inheritParams joint_sem
#' @param primary Name of the primary endpoint
#' @param n_boot Number of bootstrap replications to perform
#' @param ci_level The confidence level required
#' @param categorical Character vector of endpoint names to be treated as
#'   categorical (binary or ordinal) variables
#' @return A data.frame with point estimates, standard error estimates, and
#'   bootstrapped percentile confidence interval bounds
#' @importFrom stats sd quantile
#' @examples
#' \donttest{
#' data(joint_example)
#' # Using n_boot = 5 for example purposes only, use more bootstrap iterations
#' # in practice!
#' joint_ma_sim_3(
#'   data0 = joint_example, endpoints = c("Y1", "Y2", "Y3"), treatment = "A",
#'   primary = "Y1", n_boot = 5)
#' }
#'
#' @export
joint_ma_sim_3 <- function(
    data0, endpoints, categorical = c(), treatment = "A", primary, n_boot = 100,
    ci_level = 0.95, sandwich = FALSE, ...) {

  # Estimates on data0
  estimate <- estimate_ma_3(data0, endpoints, categorical, treatment, primary)

  # Bootstrapping
  n <- nrow(data0)
  boot_results <- as.data.frame(t(replicate(
    n_boot,
    expr = {
      data_boot <- data0[sample(1:n, size = n, replace = TRUE), ]
      estimate_ma_3(data_boot, endpoints, categorical, treatment, primary)
    }
  )))

  # Results
  re <- data.frame(
    Method = c("SEM", "Saturated", "SL", "BIC"),
    Estimate = estimate[c("est_sem", "est_saturated", "est_sl", "est_bic")],
    se = sqrt(c(estimate[c("v_sem", "v_saturated")], NA, NA)),
    se_boot = c(sd(boot_results$est_sem), sd(boot_results$est_saturated),
                sd(boot_results$est_sl), sd(boot_results$est_bic)),
    lb_boot = c(
      quantile(boot_results$est_sem, (1 - ci_level) / 2),
      quantile(boot_results$est_saturated, (1 - ci_level) / 2),
      quantile(boot_results$est_sl, (1 - ci_level) / 2),
      quantile(boot_results$est_bic, (1 - ci_level) / 2)
    ),
    ub_boot = c(
      quantile(boot_results$est_sem, (1 + ci_level) / 2),
      quantile(boot_results$est_saturated, (1 + ci_level) / 2),
      quantile(boot_results$est_sl, (1 + ci_level) / 2),
      quantile(boot_results$est_bic, (1 + ci_level) / 2)
    ),
    Omega = c(NA, NA, estimate[c("omega_sl", "omega_bic")])
  )
  rownames(re) <- NULL

  return(re)
}

#' Estimate the ATE with various methods for simulation studies
#'
#' @inheritParams joint_sem
#' @param primary Name of the primary endpoint
#' @param n_boot Number of bootstrap replications to perform
#' @param ci_level The confidence level required
#' @param categorical Character vector of endpoint names to be treated as
#'   categorical (binary or ordinal) variables
#' @return A data.frame with point estimates, standard error estimates, and
#'   bootstrapped percentile confidence interval bounds
#' @import sl3
#' @examples
#' \donttest{
#' data(joint_example)
#' # Using n_boot = 5 for example purposes only, use more bootstrap iterations
#' # in practice!
#' joint_ma_sim_4(
#'   data0 = joint_example, endpoints = c("Y1", "Y2", "Y3_cat"), treatment = "A",
#'   primary = "Y3_cat", n_boot = 5)
#' }
#' @export
joint_ma_sim_4 <- function(
    data0, endpoints, categorical = c(), treatment = "A", primary, n_boot = 100,
    ci_level = 0.95, sandwich = FALSE, ...) {

  # Estimates on data0
  estimate <- estimate_ma_4(data0, endpoints, categorical, treatment, primary)

  # Bootstrapping
  n <- nrow(data0)
  boot_results <- as.data.frame(t(replicate(
    n_boot,
    expr = {
      data_boot <- data0[sample(1:n, size = n, replace = TRUE), ]

      tryCatch(
        expr = {
          estimate_ma_4(data_boot, endpoints, categorical, treatment, primary)
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
      max(boot_results$nlm_iterations_sem, na.rm = TRUE),
      max(boot_results$nlm_iterations_saturated, na.rm = TRUE),
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


#' Estimate the ATE with various methods for simulation studies
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
#' # Continuoius endpoint (ATE)
#' joint_ma_sim_5(
#'   data0 = joint_example,
#'   endpoints = c("Y1", "Y2", "Y3_cat"),
#'   treatment = "A",
#'   primary = "Y1",
#'   n_boot = 1)
#'
#' # Binary endpoint (risk difference ATE)
#' joint_ma_sim_5(
#'   data0 = joint_example,
#'   endpoints = c("Y1", "Y2", "Y3_cat"),
#'   treatment = "A",
#'   primary = "Y3_cat",
#'   n_boot = 1)
#'
#' # Ordinal endpoint (probit regression coefficient)
#' joint_ma_sim_5(
#'   data0 = joint_example,
#'   endpoints = c("Y1", "Y2", "Y4_cat"),
#'   treatment = "A",
#'   primary = "Y4_cat",
#'   n_boot = 1)
#'
#' }
#' @export
joint_ma_sim_5 <- function(
    data0, endpoints, treatment = "A", primary, n_boot = 100,
    ci_level = 0.95, sandwich = FALSE, ...) {


  # Estimates on data0
  estimate <- estimate_ma_5(data0, endpoints, treatment, primary, ...)

  # Bootstrapping
  n <- nrow(data0)
  boot_results <- as.data.frame(t(replicate(
    n_boot,
    expr = {
      data_boot <- data0[sample(1:n, size = n, replace = TRUE), ]

      tryCatch(
        expr = {
          estimate_ma_5(data_boot, endpoints, treatment, primary)
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

sem_converge_sim <- function(
    data0, endpoints, categorical = c(), treatment = "A",
    phi_init = NULL, ...) {

}
