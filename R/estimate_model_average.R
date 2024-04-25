## FUNCTIONS FOR MODEL AVERAGING

#' Estimate the ATE with various methods for simulation studies
#'
#' @inheritParams joint_sem
#' @return A data.frame with point estimates, standard error estimates, test
#'   statistics, and 95% confidence interval bounds
#' @export
joint_ma_sim <- function(
    data0, endpoints, categorical = c(), treatment = "A", nboot = 100, ...) {


  est_ma <- estimate_ma(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, ...)

  # Estimate sampling distribution of ma estimator treating weights as fixed
  v_ma_fixed <- drop(omega %*% (est_ma$ate_ma^2 + est_ma$v_models)) -
    est_ma$ate_ma^2

  Z_ma_fixed <- est_ma$ate_ma / sqrt(v_ma_fixed)

  alpha <- 0.05

  # Estimated CDF
  F_ma <- function(x) {
    drop(est_ma$omega %*%
      pnorm(x, mean = est_ma$ate_models, sd = sqrt(est_ma$v_models))
    )
  }

  # Estimate %tile CI from CDF
  lb_fixed <- uniroot(
    f = function(x) F_ma(x) - alpha / 2,
    interval = est_ma$ate_ma + c(-5, 5) * sqrt(v_ma)
    )$root
  ub_fixed <- uniroot(
    f = function(x) F_ma(x) - (1 - alpha / 2),
    interval = est_ma$ate_ma + c(-5, 5) * sqrt(v_ma)
  )$root


  # Estimate sampling distribution via bootstrap

  boot_ests <- replicate(
    n_boot, do_boot_ma(data0, endpoints, categorical, treatment, ...)$ate_ma
    )
  v_ma_boot <- var(boot_ests)
  Z_ma_boot <- est_ma$ate_ma / sqrt(v_ma_boot)
  lb_boot <- quantile(boot_ests, alpha / 2)
  ub_boot <- quantile(boot_ests, 1 - alpha / 2)


  re <- data.frame(
    Method = c("SEM", "Saturated", "MA Fixed", "MA Boot"),
    estimate = c(est_ma$ate_models, rep(est_ma$ate_ma, 2)),
    se = sqrt(c(est_ma$v_models, v_ma_fixed, v_ma_boot))
  )
  re$Z <- re$estimate / re$se
  re$lb <- c(
    qnorm(alpha / 2, mean = est_ma$ate_models, sd = sqrt(est_ma$v_models)),
    lb_fixed, lb_boot
    )
  re$ub <- c(
    qnorm(1 - alpha / 2, mean = est_ma$ate_models, sd = sqrt(est_ma$v_models)),
    ub_fixed, ub_boot
  )
  return(re)
}


do_boot_ma <- function(data0, endpoints, categorical, treatment, ...) {

  n <- nrow(data0)
  data_boot <- data0[sample(1:n, size = n, replace = TRUE), ]

  estimate_ma(data_boot, endpoints, categorical, treatment, ...)
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
