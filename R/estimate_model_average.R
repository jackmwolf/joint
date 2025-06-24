## FUNCTIONS FOR MODEL AVERAGING ===============================================

#' Estimate the ATE and conditional means using the saturated, SEM, BMA, and SL
#' estimators
#'
#' @inheritParams joint_sem
#' @param primary Name of the primary endpoint
#' @param n_boot Number of bootstrap samples
#' @param ci_level Confidence level
#' @param force_probit For binary endpoints only. If TRUE, reports probit
#'   regression coefficient instead of risk difference.
#' @importFrom stats lm coef vcov sd quantile
#' @importFrom rpart rpart.control
#' @importFrom VGAM vglm cumulative
#' @return A vector of point estimates using the saturated, SEM, BMA, and SL
#'   estimators, variance estimates, and BMA and SL model weights
estimate_ma <- function(
    data0, endpoints, treatment = "A", primary, n_boot = 100,
    ci_level = 0.95, sandwich = FALSE, force_probit = FALSE, ...) {

  # List of categorical (factor) endpoints
  categorical <- endpoints[which(sapply(data0[, endpoints], is.factor))]

  # Temporarily code secondary categorical endpoints as numeric to avoid sl3
  # coding as indicators.
  data0_sl <- data0
  for (p in setdiff(categorical, primary)) {
    data0_sl[[p]] <- as.numeric(data0[[p]])
  }

  primary_ordinal <- is.factor(data0_sl[[primary]]) &
    length(levels(data0_sl[[primary]])) > 2

  # Initialize super learning
  secondary <- endpoints[endpoints != primary]
  task <- sl3::make_sl3_Task(
    data = data0_sl,
    outcome = primary,
    covariates = c(secondary, treatment)
  )

  # SEM learner
  lrnr_sem <- Lrnr_sem$new(
    name = "SEM", categorical = categorical, treatment = treatment, ...)

  # Define appropriate saturated learner
  if (primary_ordinal) {
    lrnr_saturated <- Lrnr_vglm$new(covariates = treatment, treatment = treatment)

  } else {
    lrnr_saturated <- sl3::Lrnr_glm$new(covariates = treatment)
  }

  stack <- sl3::Stack$new(lrnr_sem, lrnr_saturated)

  # Fit Super Learner and sub models
  if (primary_ordinal) {

    sl <- sl3::Lrnr_sl$new(learners = stack)

  } else {

    metalearner <- sl3::make_learner(
      sl3::Lrnr_solnp,
      eval_function = sl3::loss_squared_error,
      learner_function = sl3::metalearner_linear
    )

    sl <- sl3::Lrnr_sl$new(learners = stack, metalearner = metalearner)

  }

  sl_fit <- sl$train(task = task)

  # Sub-models ---
  # Extract joint_sem() model from sl_fit
  est_sem <- sl_fit$learner_fits[["SEM"]]$fit_object
  xyz_sem <- estimate_effects(est_sem, sandwich = sandwich, risk_difference = categorical)
  mu_sem <- estimate_means(est_sem, sandwich = sandwich)


  # Fit saturated model
  est_saturated <- joint_saturated(
    data0 = data0, endpoints = endpoints, treatment = treatment
  )

  # Pull estimated effects
  if (primary_ordinal | force_probit) {
    ate_sem <- unname(xyz_sem$estimate[paste0(primary, "_SD")])
    v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary, "_SD")])

    xyz_saturated <- VGAM::vglm(
      as.numeric(data0[[primary]]) ~ data0[[treatment]],
      family = VGAM::cumulative(link = "probitlink", parallel = TRUE)
    )

    ate_saturated <- unname(-1 * coef(xyz_saturated)[length(coef(xyz_saturated))])
    v_saturated <- unname(diag(vcov(xyz_saturated))[length(coef(xyz_saturated))])

    mu0_satuated <- NA
    mu1_satuated <- NA

    v_mu0_saturated <- NA
    v_mu1_saturated <- NA

    mu0_sem <- NA
    mu1_sem <- NA

    v_mu0_sem <- NA
    v_mu1_sem <- NA
    cov_mu0_mu1_sem <- NA

  } else {

    if (primary %in% categorical) {
      ate_sem <- unname(xyz_sem$estimate[paste0(primary, "_RD")])
      v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary, "_RD")])
    } else {
      ate_sem <- unname(xyz_sem$estimate[paste0(primary)])
      v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary)])
    }

    xyz_saturated <- lm(as.numeric(data0[[primary]]) ~ data0[[treatment]])
    ate_saturated <- unname(coef(xyz_saturated)[2])
    v_saturated <- unname(vcov(xyz_saturated)[2, 2])

    mu0_satuated <- mean(as.numeric(data0[[primary]][data0[[treatment]] == 0]))
    mu1_satuated <- mean(as.numeric(data0[[primary]][data0[[treatment]] == 1]))

    v_mu0_saturated <- var(as.numeric(data0[[primary]][data0[[treatment]] == 0])) /
      sum(data0[[treatment]] == 0)
    v_mu1_saturated <- var(as.numeric(data0[[primary]][data0[[treatment]] == 1])) /
      sum(data0[[treatment]] == 1)

    mu0_sem <- unname(mu_sem$estimate[paste0(primary, "_A0")])
    mu1_sem <- unname(mu_sem$estimate[paste0(primary, "_A1")])

    v_mu0_sem <- unname(diag(mu_sem$vcov)[paste0(primary, "_A0")])
    v_mu1_sem <- unname(diag(mu_sem$vcov)[paste0(primary, "_A1")])
    cov_mu0_mu1_sem <-
      unname(mu_sem$vcov[paste0(primary, "_A0"), paste0(primary, "_A1")])

  }

  # BIC Weights ---
  # log likelihoods
  lls <- c(est_sem$ll, est_saturated$ll)
  names(lls) <- c("sem", "saturated")

  # Number of parameters
  pk <- c(est_sem$dim_phi, est_saturated$dim_phi)

  # BICs
  bics <- -2 * lls + pk * log(nrow(data0))

  # Weights
  omega_bic <-
    exp(-1/2 * (bics - min(bics))) / sum(exp(-1/2 * (bics - min(bics))))
  ate_ma_bic <- drop(omega_bic %*% c(ate_sem, ate_saturated))

  # SL Weights ---
  omega_sl <- sl_fit$coefficients
  ate_ma_sl <- drop(omega_sl %*% c(ate_sem, ate_saturated))

  # Return point estimates, sub-model variance estimates, and weights ---
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
    nlm_iterations_saturated = est_saturated$nlm_iterations,

    mu0_satuated = mu0_satuated,
    mu1_satuated = mu1_satuated,
    v_mu0_saturated = v_mu0_saturated,
    v_mu1_saturated = v_mu1_saturated,
    mu0_sem = mu0_sem,
    mu1_sem = mu1_sem,
    v_mu0_sem = v_mu0_sem,
    v_mu1_sem = v_mu1_sem,
    cov_mu0_mu1_sem = cov_mu0_mu1_sem

  )

}

#' Estimate the ATE using the SL-SEM estimator
#'
#' \code{estimate_sl} estimates the ATE and conditional means of the primary
#' endpoint on one dataset using super learning.
#'
#' @inheritParams joint_sl
#'
#' @importFrom stats lm coef vcov
#' @importFrom rpart rpart.control
#' @importFrom VGAM vglm cumulative
#'
#' @return A vector of point estimates using the saturated, SEM, and SL
#'   estimators, variance estimates, and and SL model weights
estimate_sl <- function(
    data0, endpoints, treatment = "A", primary,
    ci_level = 0.95, force_probit = FALSE, ...) {
  # List of categorical (factor) endpoints
  categorical <- endpoints[which(sapply(data0[, endpoints], is.factor))]

  # Temporarily code secondary categorical endpoints as numeric to avoid sl3
  # coding as indicators.
  data0_sl <- data0
  for (p in setdiff(categorical, primary)) {
    data0_sl[[p]] <- as.numeric(data0[[p]])
  }

  primary_ordinal <- is.factor(data0_sl[[primary]]) &
    length(levels(data0_sl[[primary]])) > 2

  # Initialize super learning
  secondary <- endpoints[endpoints != primary]
  task <- sl3::make_sl3_Task(
    data = data0_sl,
    outcome = primary,
    covariates = c(secondary, treatment)
  )

  # SEM learner
  lrnr_sem <- Lrnr_sem$new(
    name = "SEM", categorical = categorical, treatment = treatment, ...)

  # Define appropriate saturated learner
  if (primary_ordinal) {
    lrnr_saturated <- Lrnr_vglm$new(covariates = treatment, treatment = treatment)

  } else {
    lrnr_saturated <- sl3::Lrnr_glm$new(covariates = treatment)
  }

  stack <- sl3::Stack$new(lrnr_sem, lrnr_saturated)

  # Fit Super Learner and sub models
  if (primary_ordinal) {

    sl <- sl3::Lrnr_sl$new(learners = stack)

  } else {

    metalearner <- sl3::make_learner(
      sl3::Lrnr_solnp,
      eval_function = sl3::loss_squared_error,
      learner_function = sl3::metalearner_linear
    )

    sl <- sl3::Lrnr_sl$new(learners = stack, metalearner = metalearner)

  }

  sl_fit <- sl$train(task = task)

  # Sub-models ---
  # Extract joint_sem() model from sl_fit
  est_sem <- sl_fit$learner_fits[["SEM"]]$fit_object
  xyz_sem <- estimate_effects(est_sem, risk_difference = categorical)
  mu_sem <- estimate_means(est_sem)


  # Fit saturated model
  est_saturated <- joint_saturated(
    data0 = data0, endpoints = endpoints, treatment = treatment
  )

  # Pull estimated effects
  if (primary_ordinal | force_probit) {
    ate_sem <- unname(xyz_sem$estimate[paste0(primary, "_SD")])
    v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary, "_SD")])

    xyz_saturated <- VGAM::vglm(
      as.numeric(data0[[primary]]) ~ data0[[treatment]],
      family = VGAM::cumulative(link = "probitlink", parallel = TRUE)
    )

    ate_saturated <- unname(-1 * coef(xyz_saturated)[length(coef(xyz_saturated))])
    v_saturated <- unname(diag(vcov(xyz_saturated))[length(coef(xyz_saturated))])

    mu0_satuated <- NA
    mu1_satuated <- NA

    v_mu0_saturated <- NA
    v_mu1_saturated <- NA

    mu0_sem <- NA
    mu1_sem <- NA

    v_mu0_sem <- NA
    v_mu1_sem <- NA
    cov_mu0_mu1_sem <- NA

  } else {

    if (primary %in% categorical) {
      ate_sem <- unname(xyz_sem$estimate[paste0(primary, "_RD")])
      v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary, "_RD")])
    } else {
      ate_sem <- unname(xyz_sem$estimate[paste0(primary)])
      v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary)])
    }

    xyz_saturated <- lm(as.numeric(data0[[primary]]) ~ data0[[treatment]])
    ate_saturated <- unname(coef(xyz_saturated)[2])
    v_saturated <- unname(vcov(xyz_saturated)[2, 2])

    mu0_satuated <- mean(as.numeric(data0[[primary]][data0[[treatment]] == 0]))
    mu1_satuated <- mean(as.numeric(data0[[primary]][data0[[treatment]] == 1]))

    v_mu0_saturated <- var(as.numeric(data0[[primary]][data0[[treatment]] == 0])) /
      sum(data0[[treatment]] == 0)
    v_mu1_saturated <- var(as.numeric(data0[[primary]][data0[[treatment]] == 1])) /
      sum(data0[[treatment]] == 1)

    mu0_sem <- unname(mu_sem$estimate[paste0(primary, "_A0")])
    mu1_sem <- unname(mu_sem$estimate[paste0(primary, "_A1")])

    v_mu0_sem <- unname(diag(mu_sem$vcov)[paste0(primary, "_A0")])
    v_mu1_sem <- unname(diag(mu_sem$vcov)[paste0(primary, "_A1")])
    cov_mu0_mu1_sem <-
      unname(mu_sem$vcov[paste0(primary, "_A0"), paste0(primary, "_A1")])

  }

  # SL Weights ---
  omega_sl <- sl_fit$coefficients
  ate_ma_sl <- drop(omega_sl %*% c(ate_sem, ate_saturated))

  # Return point estimates, sub-model variance estimates, and weights ---
  c(
    est_sl        = ate_ma_sl,
    est_sem       = ate_sem,
    est_saturated = ate_saturated,

    omega_sl  = unname(omega_sl[1]),

    v_sem       = v_sem,
    v_saturated = v_saturated,

    nlm_code_sem = est_sem$nlm_code,
    nlm_iterations_sem = est_sem$nlm_iterations,

    nlm_code_saturated = est_saturated$nlm_code,
    nlm_iterations_saturated = est_saturated$nlm_iterations,

    mu0_satuated = mu0_satuated,
    mu1_satuated = mu1_satuated,
    v_mu0_saturated = v_mu0_saturated,
    v_mu1_saturated = v_mu1_saturated,
    mu0_sem = mu0_sem,
    mu1_sem = mu1_sem,
    v_mu0_sem = v_mu0_sem,
    v_mu1_sem = v_mu1_sem,
    cov_mu0_mu1_sem = cov_mu0_mu1_sem
  )
}

#' Estimate the ATE using the BMA-SEM estimator
#'
#' \code{estimate_bma} estimates the ATE and conditional means of the primary
#' endpoint on one dataset using Bayesian model averaging.
#'
#' @inheritParams joint_bma
#' @inherit estimate_ma return
#'
#' @importFrom stats lm coef vcov
#' @importFrom rpart rpart.control
#' @importFrom VGAM vglm cumulative
#'
#' @return A vector of point estimates using the saturated, SEM, and BMA
#'   estimators, variance estimates, and and BMA model weights
estimate_bma <- function(
    data0, endpoints, treatment = "A", primary,
    ci_level = 0.95, force_probit = FALSE, ...) {

  # List of categorical (factor) endpoints
  categorical <- endpoints[which(sapply(data0[, endpoints], is.factor))]

  primary_ordinal <- is.factor(data0[[primary]]) &
    length(levels(data0[[primary]])) > 2

  est_sem <- joint_sem(
    data0 = data0,
    endpoints = endpoints,
    treatment = treatment,
    ...
  )
  est_saturated <- joint_saturated(
    data0 = data0,
    endpoints = endpoints,
    treatment = treatment
  )

  xyz_sem <- estimate_effects(est_sem, sandwich = FALSE, risk_difference = categorical)
  mu_sem <- estimate_means(est_sem, sandwich = FALSE)

  # Pull estimated effects
  if (primary_ordinal | force_probit) {
    ate_sem <- unname(xyz_sem$estimate[paste0(primary, "_SD")])
    v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary, "_SD")])

    xyz_saturated <- VGAM::vglm(
      as.numeric(data0[[primary]]) ~ data0[[treatment]],
      family = VGAM::cumulative(link = "probitlink", parallel = TRUE)
    )

    ate_saturated <- unname(-1 * coef(xyz_saturated)[length(coef(xyz_saturated))])
    v_saturated <- unname(diag(vcov(xyz_saturated))[length(coef(xyz_saturated))])

    mu0_satuated <- NA
    mu1_satuated <- NA

    v_mu0_saturated <- NA
    v_mu1_saturated <- NA

    mu0_sem <- NA
    mu1_sem <- NA

    v_mu0_sem <- NA
    v_mu1_sem <- NA
    cov_mu0_mu1_sem <- NA

  } else {

    if (primary %in% categorical) {
      ate_sem <- unname(xyz_sem$estimate[paste0(primary, "_RD")])
      v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary, "_RD")])
    } else {
      ate_sem <- unname(xyz_sem$estimate[paste0(primary)])
      v_sem  <- unname(diag(xyz_sem$vcov)[paste0(primary)])
    }

    xyz_saturated <- lm(as.numeric(data0[[primary]]) ~ data0[[treatment]])
    ate_saturated <- unname(coef(xyz_saturated)[2])
    v_saturated <- unname(vcov(xyz_saturated)[2, 2])

    mu0_satuated <- mean(as.numeric(data0[[primary]][data0[[treatment]] == 0]))
    mu1_satuated <- mean(as.numeric(data0[[primary]][data0[[treatment]] == 1]))

    v_mu0_saturated <- var(as.numeric(data0[[primary]][data0[[treatment]] == 0])) /
      sum(data0[[treatment]] == 0)
    v_mu1_saturated <- var(as.numeric(data0[[primary]][data0[[treatment]] == 1])) /
      sum(data0[[treatment]] == 1)

    mu0_sem <- unname(mu_sem$estimate[paste0(primary, "_A0")])
    mu1_sem <- unname(mu_sem$estimate[paste0(primary, "_A1")])

    v_mu0_sem <- unname(diag(mu_sem$vcov)[paste0(primary, "_A0")])
    v_mu1_sem <- unname(diag(mu_sem$vcov)[paste0(primary, "_A1")])
    cov_mu0_mu1_sem <-
      unname(mu_sem$vcov[paste0(primary, "_A0"), paste0(primary, "_A1")])

  }


  # BIC Weights ---
  # log likelihoods
  lls <- c(est_sem$ll, est_saturated$ll)
  names(lls) <- c("sem", "saturated")

  # Number of parameters
  pk <- c(est_sem$dim_phi, est_saturated$dim_phi)

  # BICs
  bics <- -2 * lls + pk * log(nrow(data0))

  # Weights
  omega_bic <-
    exp(-1/2 * (bics - min(bics))) / sum(exp(-1/2 * (bics - min(bics))))
  ate_ma_bic <- drop(omega_bic %*% c(ate_sem, ate_saturated))

  # Return point estimates, sub-model variance estimates, and weights ---
  c(
    est_bic       = ate_ma_bic,
    est_sem       = ate_sem,
    est_saturated = ate_saturated,

    omega_bic = unname(omega_bic[1]),

    v_sem       = v_sem,
    v_saturated = v_saturated,

    nlm_code_sem = est_sem$nlm_code,
    nlm_iterations_sem = est_sem$nlm_iterations,

    nlm_code_saturated = est_saturated$nlm_code,
    nlm_iterations_saturated = est_saturated$nlm_iterations,

    mu0_satuated = mu0_satuated,
    mu1_satuated = mu1_satuated,
    v_mu0_saturated = v_mu0_saturated,
    v_mu1_saturated = v_mu1_saturated,
    mu0_sem = mu0_sem,
    mu1_sem = mu1_sem,
    v_mu0_sem = v_mu0_sem,
    v_mu1_sem = v_mu1_sem,
    cov_mu0_mu1_sem = cov_mu0_mu1_sem

  )
}

#' Estimate the ATE using the SL-SEM estimator
#'
#' \code{joint_sl} estimates the ATE and conditional means of the primary
#' endpoint by super learning between the saturated model and one-factor
#' structural equation model.
#'
#' @inheritParams estimate_ma
#'
#' @return A list containing the following named objects:
#'     \item{summary}{A table with summarized output including ATE estimates,
#'     standard error estimates, and bootstrapped percentile CI bounds}
#'     \item{estimate}{A vector of parameter estimates and diagnostic
#'     information}
#'     \item{boot_results}{A \code{data.frame} where each row is an estimate
#'     obtained from a particular bootstrapped resample of \code{data0}}
#'     \item{runtime}{The run time in seconds}
#'     \item{call}{The matched call}
#'
#' @importFrom Rdpack reprompt
#' @references{ \insertRef{wolf_jointly_2025}{joint} }
#'
#' @seealso [joint_bma()]
#'
#' @examples
#' \donttest{
#' data(joint_example)
#' set.seed(1)
#' joint_sl(
#'   data0 = joint_example,
#'   endpoints = c("Y1", "Y2", "Y3"),
#'   treatment = "A",
#'   primary = "Y1",
#'   n_boot = 5,
#'   ci_level = 0.95
#'   )
#' }
#'
#' @export
#'
joint_sl <- function(
    data0, endpoints, treatment = "A", primary, n_boot = 100,
    ci_level = 0.95, force_probit = FALSE, ...) {

  cl <- match.call()
  t0 <- Sys.time()

  # Verify input
  verify_input_joint_sl(
    data0, endpoints, treatment, primary, n_boot, ci_level
  )

  # Estimates on data0
  estimate <- estimate_sl(data0, endpoints, treatment, primary, ...)

  # Bootstrapping
  n <- nrow(data0)
  boot_results <- as.data.frame(t(replicate(
    n_boot,
    expr = {
      data_boot <- data0[sample(1:n, size = n, replace = TRUE), ]

      tryCatch(
        expr = {
          estimate_sl(data_boot, endpoints, treatment, primary, ...)
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
  re_summary <- data.frame(
    Method = c("SL-SEM", "SEM", "Saturated"),

    Estimate = estimate[c("est_sl", "est_sem", "est_saturated")],

    se = sqrt(c(NA, estimate[c("v_sem", "v_saturated")])),

    se_boot = c(
      sd(boot_results$est_sl, na.rm = TRUE),
      sd(boot_results$est_sem, na.rm = TRUE),
      sd(boot_results$est_saturated, na.rm = TRUE)
    ),

    lb_boot = c(
      quantile(boot_results$est_sl, (1 - ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_sem, (1 - ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_saturated, (1 - ci_level) / 2, na.rm = TRUE)
    ),

    ub_boot = c(
      quantile(boot_results$est_sl, (1 + ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_sem, (1 + ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_saturated, (1 + ci_level) / 2, na.rm = TRUE)
    ),

    Omega = c(estimate[c("omega_sl")], 1, 0)

  )
  rownames(re_summary) <- NULL

  t1 <- Sys.time()

  re <- list(
    summary = re_summary,
    estimate = estimate,
    boot_results = boot_results,
    runtime = as.numeric(t1 - t0, units = "secs")
  )
  re$call <- cl

  return(re)
}

#' Estimate the ATE using the BMA-SEM estimator
#'
#' \code{joint_bma} estimates the ATE and conditional means of the primary
#' endpoint by Bayesian model averaging between the saturated model and
#' one-factor structural equation model.
#'
#' @inheritParams joint_sl
#' @inherit joint_sl return
#'
#' @importFrom Rdpack reprompt
#' @references{ \insertRef{wolf_jointly_2025}{joint} }
#'
#' @seealso [joint_sl()]
#'
#' @examples
#' \donttest{
#' data(joint_example)
#' set.seed(1)
#' joint_bma(
#'   data0 = joint_example,
#'   endpoints = c("Y1", "Y2", "Y3"),
#'   treatment = "A",
#'   primary = "Y1",
#'   n_boot = 5,
#'   ci_level = 0.95
#'   )
#' }
#'
#' @export
joint_bma <- function(
    data0, endpoints, treatment = "A", primary, n_boot = 100,
    ci_level = 0.95, force_probit = FALSE, ...) {

  cl <- match.call()
  t0 <- Sys.time()

  # Verify input
  verify_input_joint_sl(
    data0, endpoints, treatment, primary, n_boot, ci_level
  )

  # Estimates on data0
  estimate <- estimate_bma(data0, endpoints, treatment, primary, ...)

  # Bootstrapping
  n <- nrow(data0)
  boot_results <- as.data.frame(t(replicate(
    n_boot,
    expr = {
      data_boot <- data0[sample(1:n, size = n, replace = TRUE), ]

      tryCatch(
        expr = {
          estimate_bma(data_boot, endpoints, treatment, primary, ...)
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
  re_summary <- data.frame(
    Method = c("SL-BMA", "SEM", "Saturated"),

    Estimate = estimate[c("est_bic", "est_sem", "est_saturated")],

    se = sqrt(c(NA, estimate[c("v_sem", "v_saturated")])),

    se_boot = c(
      sd(boot_results$est_bic, na.rm = TRUE),
      sd(boot_results$est_sem, na.rm = TRUE),
      sd(boot_results$est_saturated, na.rm = TRUE)
    ),

    lb_boot = c(
      quantile(boot_results$est_bic, (1 - ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_sem, (1 - ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_saturated, (1 - ci_level) / 2, na.rm = TRUE)
    ),

    ub_boot = c(
      quantile(boot_results$est_bic, (1 + ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_sem, (1 + ci_level) / 2, na.rm = TRUE),
      quantile(boot_results$est_saturated, (1 + ci_level) / 2, na.rm = TRUE)
    ),

    Omega = c(estimate[c("omega_bic")], 1, 0)

  )
  rownames(re_summary) <- NULL

  t1 <- Sys.time()

  re <- list(
    summary = re_summary,
    estimate = estimate,
    boot_results = boot_results,
    runtime = as.numeric(t1 - t0, units = "secs")
  )
  re$call <- cl

  return(re)
}
