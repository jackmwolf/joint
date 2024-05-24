## FUNCTIONS FOR MODEL AVERAGING

#' Estimate the ATE with various methods for simulation studies
#'
#' @inheritParams joint_sem
#' @param n_boot Number of bootstrap replications to perform
#' @param ci_level The confidence level required
#' @return A data.frame with point estimates, standard error estimates, and
#'   bootstrapped percentile confidence interval bounds
#' @examples
#' data(joint_example)
#' # Using n_boot = 5 for example purposes only, use more bootstrap iterations
#' # in practice!
#' joint_ma_sim_3(
#'   data0 = joint_example, endpoints = c("Y1", "Y2", "Y3"), treatment = "A",
#'   primary = "Y1", n_boot = 5)
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
#' @param n_boot Number of bootstrap replications to perform
#' @param ci_level The confidence level required
#' @return A data.frame with point estimates, standard error estimates, and
#'   bootstrapped percentile confidence interval bounds
#' @import sl3
#' @examples
#' data(joint_example)
#' # Using n_boot = 5 for example purposes only, use more bootstrap iterations
#' # in practice!
#' joint_ma_sim_4(
#'   data0 = joint_example, endpoints = c("Y1", "Y2", "Y3_cat"), treatment = "A",
#'   primary = "Y3_cat", n_boot = 0)
#'
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
      estimate_ma_4(data_boot, endpoints, categorical, treatment, primary)
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
    Omega = c(NA, NA, estimate[c("omega_sl", "omega_bic")]),
    iterations = c(
      estimate[c("nlm_iterations_sem", "nlm_iterations_saturated")],
      NA,
      NA),
    max_iterations_boot = c(
      max(boot_results$nlm_iterations_sem),
      max(boot_results$nlm_iterations_saturated),
      NA,
      NA
    ),
    nlm_code =c(
      estimate[c("nlm_code_sem", "nlm_code_saturated")],
      NA,
      NA)
  )
  rownames(re) <- NULL

  return(re)
}

estimate_ma_3 <- function(
    data0, endpoints, categorical = c(), treatment = "A", primary, n_boot = 100,
    ci_level = 0.95, sandwich = FALSE, ...) {

  # Initialize super learning
  secondary <- endpoints[endpoints != primary]
  task <- sl3::make_sl3_Task(
    data = data0,
    outcome = primary,
    covariates = c(secondary, "A")
  )

  lrnr_sem <- Lrnr_sem$new(categorical = categorical, ...)
  lrnr_saturated <- sl3::Lrnr_glm$new(covariates = treatment)

  stack <- sl3::Stack$new(lrnr_sem, lrnr_saturated)

  # Fit Super Learner and sub models
  sl <- sl3::Lrnr_sl$new(learners = stack)
  sl_fit <- sl$train(task = task)

  # Sub-models ---
  # Extract joint_sem() model from sl_fit
  est_sem <- sl_fit$learner_fits$Lrnr_sem_NULL$fit_object
  xyz_sem <- estimate_effects(est_sem, sandwich = sandwich)
  ate_sem <- unname(xyz_sem$estimate[primary])
  v_sem  <- unname(diag(xyz_sem$vcov)[primary])

  # Fit saturated model
  est_saturated <- joint_saturated(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, ...
  )
  xyz_saturated <- lm(data0[[primary]] ~ data0[[treatment]])
  ate_saturated <- unname(coef(xyz_saturated)[2])
  v_saturated <- unname(vcov(xyz_saturated)[2, 2])

  # BIC Weights ---
  # log likelihoods
  lls <- c(est_sem$ll, est_saturated$ll)
  names(lls) <- c("sem", "saturated")

  # Number of parameters
  pk <- c(est_sem$dim_phi, est_saturated$dim_phi)

  # BICs
  bics <- -2 * lls + pk * log(nrow(data0))

  # Weights
  omega_bic <- exp(-1/2 * (bics - min(bics))) / sum(exp(-1/2 * (bics - min(bics))))
  ate_ma_bic <- drop(omega_bic %*% c(ate_sem, ate_saturated))
  omega_sl <- sl_fit$coefficients
  ate_ma_sl <- drop(omega_sl %*% c(ate_sem, ate_saturated))

  # Return point estimates, sub-model variance estimates, and weights
  c(
    est_sl        = ate_ma_sl,
    est_bic       = ate_ma_bic,
    est_sem       = ate_sem,
    est_saturated = ate_saturated,

    omega_sl  = unname(omega_sl[1]),
    omega_bic = unname(omega_bic[1]),

    v_sem       = v_sem,
    v_saturated = v_saturated
  )

}

estimate_ma_4 <- function(
    data0, endpoints, categorical = c(), treatment = "A", primary, n_boot = 100,
    ci_level = 0.95, sandwich = FALSE, ...) {

  # Initialize super learning
  secondary <- endpoints[endpoints != primary]
  task <- sl3::make_sl3_Task(
    data = data0,
    outcome = primary,
    covariates = c(secondary, "A")
  )

  lrnr_sem <- Lrnr_sem$new(categorical = categorical, ...)
  lrnr_saturated <- sl3::Lrnr_glm$new(covariates = treatment)

  stack <- sl3::Stack$new(lrnr_sem, lrnr_saturated)

  # Fit Super Learner and sub models
  metalearner <- sl3::make_learner(
    sl3::Lrnr_solnp,
    eval_function = sl3::loss_squared_error,
    learner_function = sl3::metalearner_linear
    )
  sl <- sl3::Lrnr_sl$new(learners = stack, metalearner = metalearner)
  sl_fit <- sl$train(task = task)

  # Sub-models ---
  # Extract joint_sem() model from sl_fit
  est_sem <- sl_fit$learner_fits[[1]]$fit_object
  xyz_sem <- estimate_effects(est_sem, sandwich = sandwich, risk_difference = categorical)

  if (primary %in% categorical) {
    ate_sem <- unname(xyz_sem$estimate[paste0(primary, "_RD")])
    v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary, "_RD")])
  } else {
    ate_sem <- unname(xyz_sem$estimate[paste0(primary)])
    v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary)])
  }


  # Fit saturated model
  est_saturated <- joint_saturated(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, ...
  )
  xyz_saturated <- lm(data0[[primary]] ~ data0[[treatment]])
  ate_saturated <- unname(coef(xyz_saturated)[2])
  v_saturated <- unname(vcov(xyz_saturated)[2, 2])

  # BIC Weights ---
  # log likelihoods
  lls <- c(est_sem$ll, est_saturated$ll)
  names(lls) <- c("sem", "saturated")

  # Number of parameters
  pk <- c(est_sem$dim_phi, est_saturated$dim_phi)

  # BICs
  bics <- -2 * lls + pk * log(nrow(data0))

  # Weights
  omega_bic <- exp(-1/2 * (bics - min(bics))) / sum(exp(-1/2 * (bics - min(bics))))
  ate_ma_bic <- drop(omega_bic %*% c(ate_sem, ate_saturated))
  omega_sl <- sl_fit$coefficients
  ate_ma_sl <- drop(omega_sl %*% c(ate_sem, ate_saturated))

  # Return point estimates, sub-model variance estimates, and weights
  c(
    est_sl        = ate_ma_sl,
    est_bic       = ate_ma_bic,
    est_sem       = ate_sem,
    est_saturated = ate_saturated,

    omega_sl  = unname(omega_sl[1]),
    omega_bic = unname(omega_bic[1]),

    v_sem       = v_sem,
    v_saturated = v_saturated,

    nlm_code_sem = est_sem$nlm_code,
    nlm_iterations_sem = est_sem$nlm_iterations,

    nlm_code_saturated = est_saturated$nlm_code,
    nlm_iterations_saturated = est_saturated$nlm_iterations
  )

}


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

#' Estimate the model averaged ATE on one bootstrap sample
#'
#' @inheritParams joint_sem
#' @inherit estimate_ma return
do_boot_ma <- function(data0, endpoints, categorical, treatment, ...) {

  n <- nrow(data0)
  data_boot <- data0[sample(1:n, size = n, replace = TRUE), ]

  estimate_ma(data_boot, endpoints, categorical, treatment, sandwich = FALSE)
}

#' Estimate the ATE using model averaging
#'
#' @inheritParams joint_sem
#' @return A list with point estimates using BIC model averagin and
#'   sub-model estimates, SE estimates, and weights,
estimate_ma <- function(data0, endpoints, categorical, treatment, sandwich,
                        ...) {
  # Estimate SEM and saturated model
  est_sem <- joint_sem(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, sandwich = sandwich, ...
  )

  est_saturated <- joint_saturated(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, ...
  )


  # Estimate effects
  xyz_sem <- estimate_effects(est_sem, sandwich = sandwich)
  ate_sem <- unname(xyz_sem$estimate[endpoints[1]])
  v_sem  <- unname(diag(xyz_sem$vcov)[endpoints[1]])

  xyz_saturated <- lm(data0[[endpoints[1]]] ~ data0[[treatment]])
  ate_saturated <- unname(coef(xyz_saturated)[2])
  v_saturated <- unname(vcov(xyz_saturated)[2, 2])

  # log likelihoods
  lls <- c(est_sem$ll, est_saturated$ll)
  names(lls) <- c("sem", "saturated")


  # Number of parameters
  pk <- c(est_sem$dim_phi, est_saturated$dim_phi)

  # BICs
  bics <- -2 * lls + pk * log(nrow(data0))

  # Weights
  omega <- exp(-1/2 * (bics - min(bics))) / sum(exp(-1/2 * (bics - min(bics))))
  ate_ma <- drop(omega %*% c(ate_sem, ate_saturated))

  list(
    ate_ma = ate_ma,
    ate_models = c(sem = ate_sem, saturated = ate_saturated),
    v_models = c(sem = v_sem, saturated = v_saturated),
    omega = omega,
    lls = lls
  )
}

#' Estimate the ATE using model averaging
#'
#' @inheritParams joint_sem
#' @param Primary Name of the primary endpoint
#' @return A list with point estimates using SL and BIC model averaging,
#'   sub-model estimates, and weights.
estimate_ma_2 <- function(data0, endpoints, categorical, treatment = "A", primary, sandwich = FALSE,
                        ...) {
  # Estimate SEM and saturated model
  est_sem <- joint_sem(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, sandwich = sandwich, ...
  )

  est_saturated <- joint_saturated(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment, ...
  )


  # Estimate effects
  xyz_sem <- estimate_effects(est_sem, sandwich = sandwich)
  ate_sem <- unname(xyz_sem$estimate[endpoints[1]])
  v_sem  <- unname(diag(xyz_sem$vcov)[endpoints[1]])

  xyz_saturated <- lm(data0[[endpoints[1]]] ~ data0[[treatment]])
  ate_saturated <- unname(coef(xyz_saturated)[2])
  v_saturated <- unname(vcov(xyz_saturated)[2, 2])

  # log likelihoods
  lls <- c(est_sem$ll, est_saturated$ll)
  names(lls) <- c("sem", "saturated")

  # Number of parameters
  pk <- c(est_sem$dim_phi, est_saturated$dim_phi)

  # BICs
  bics <- -2 * lls + pk * log(nrow(data0))

  # Super Learning
  secondary <- endpoints[endpoints != primary]
  task <- sl3::make_sl3_Task(
    data = data0,
    outcome = primary,
    covariates = c(secondary, "A")
  )

  lrnr_sem <- Lrnr_sem$new(categorical = categorical, ...)
  lrnr_saturated <- sl3::Lrnr_glm$new(covariates = treatment)

  stack <- sl3::Stack$new(lrnr_sem, lrnr_saturated)

  sl <- sl3::Lrnr_sl$new(learners = stack)
  sl_fit <- sl$train(task = task)

  # Weights
  omega_bic <- exp(-1/2 * (bics - min(bics))) / sum(exp(-1/2 * (bics - min(bics))))
  ate_ma_bic <- drop(omega_bic %*% c(ate_sem, ate_saturated))
  omega_sl <- sl_fit$coefficients
  ate_ma_sl <- drop(omega_sl %*% c(ate_sem, ate_saturated))

  list(
    ate_ma_bic = ate_ma_bic,
    ate_ma_sl = ate_ma_sl,
    ate_models = c(sem = ate_sem, saturated = ate_saturated),
    v_models = c(sem = v_sem, saturated = v_saturated),
    omega_bic = omega_bic,
    omega_sl = omega_sl,
    lls = lls
  )
}

#' Estimate the ATE with various methods for simulation studies
#'
#' @inheritParams joint_sem
#' @inheritParams estimate_ma_2
#' @return A data.frame with point estimates for various methods
#' @examples
#' data(joint_example)
#' estimate_ma_2(
#'     joint_example,
#'     endpoints = c("Y1", "Y2", "Y3"),
#'     categorical = c(),
#'     treatment = "A",
#'     primary = "Y1"
#' )
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
