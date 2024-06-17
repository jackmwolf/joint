## FUNCTIONS FOR MODEL AVERAGING

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

