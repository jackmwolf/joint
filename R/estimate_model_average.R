## FUNCTIONS FOR MODEL AVERAGING

#' Estimate the ATE with various methods for simulation studies
#'
#' @inheritParams joint_sem
#' @param n_boot Number of bootstrap replications to perform
#' @param ci_level The confidence level required
#' @return A data.frame with point estimates, standard error estimates, test
#'   statistics, and 95% confidence interval bounds
#' @examples
#' data(joint_example)
#' # Using n_boot = 5 for example purposes only, use more bootstrap iterations
#' # in practice!
#' joint_ma_sim(joint_example, c("Y1", "Y2", "Y3"), n_boot = 5)
#'
#' @export
joint_ma_sim <- function(
    data0, endpoints, categorical = c(), treatment = "A", n_boot = 100,
    ci_level = 0.95, ...) {


  est_ma <- estimate_ma(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, ...)

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
  lb_fixed <- uniroot(
    f = function(x) F_ma(x) - alpha / 2,
    interval = est_ma$ate_ma + c(-5, 5) * sqrt(v_ma_fixed),
    extendInt = "yes"
    )$root
  ub_fixed <- uniroot(
    f = function(x) F_ma(x) - (1 - alpha / 2),
    interval = est_ma$ate_ma + c(-5, 5) * sqrt(v_ma_fixed),
    extendInt = "yes"
  )$root


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
  return(re)
}


do_boot_ma <- function(data0, endpoints, categorical, treatment, ...) {

  n <- nrow(data0)
  data_boot <- data0[sample(1:n, size = n, replace = TRUE), ]

  estimate_ma(data_boot, endpoints, categorical, treatment, sandwich = FALSE)
}

estimate_ma <- function(data0, endpoints, categorical, treatment,
                        ...) {
  # Estimate SEM and saturated model
  est_sem <- joint_sem(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, ...
  )

  est_saturated <- joint_saturated(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, ...
  )

  # Estimate effects
  xyz_sem <- estimate_effects(est_sem)
  ate_sem <- unname(xyz_sem$estimate[endpoints[1]])
  v_sem  <- unname(diag(xyz_sem$vcov)[endpoints[1]])

  xyz_saturated <- lm(data0[[endpoints[1]]] ~ data0[[treatment]])
  ate_saturated <- unname(coef(xyz_saturated)[2])
  v_saturated <- unname(vcov(xyz_saturated)[2, 2])

  # log likelihoods
  lls <- c(est_sem$ll, est_saturated$ll)

  # Number of parameters
  pk <- c(length(est_sem$estimate), est_saturated$dim_phi)

  # BICs
  bics <- -2 * lls + pk * log(nrow(data0))

  # Weights
  omega <- exp(-1/2 * (bics - min(bics)))/sum(exp(-1/2 * (bics - min(bics))))
  ate_ma <- drop(omega %*% c(ate_sem, ate_saturated))

  list(
    ate_ma = ate_ma,
    ate_models = c(sem = ate_sem, saturated = ate_saturated),
    v_models = c(sem = v_sem, saturated = v_saturated),
    omega = omega
  )
}
