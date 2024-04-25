

#' Estimate an average treatment effect using model averaging (sl3)
#'
#' @inheritParams joint_sem
#' @param ... Additional arguments to \code{joint_sem}
#' @export
#' @examples
#' data(joint_example)
#' estimate_ma <- estimate_sem_ma(
#'   data0 = joint_example,
#'   endpoints = c("Y1", "Y2", "Y3"),
#'   categorical = c(),
#'   treatment = "A",
#'   primary = "Y1"
#' )
#' estimate_ma$tmle_est
#' estimate_ma$tmle_sigma
#' estimate_sem <- estimate_effects(estimate_ma$joint_sem)
#' estimate_sem$estimate["Y1"]
#' sqrt(diag(estimate_sem$vcov)["Y1"])
estimate_sem_ma <- function(data0, endpoints, categorical = c(), treatment = "A", primary, ...) {

  secondary <- endpoints[endpoints != primary]

  # Set up Super Learning
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
  sl_fit

  # Tasks for counterfactual prediction
  dataA1 <- dataA0 <- data0
  dataA1[[treatment]] <- 1
  dataA0[[treatment]] <- 0

  task_1 <- sl3::make_sl3_Task(
    data = dataA1,
    outcome = primary,
    covariates = c(secondary, treatment)
  )
  task_0 <- sl3::make_sl3_Task(
    data = dataA0,
    outcome = primary,
    covariates = c(secondary, treatment)
  )

  Q0_1 <- sl_fit$predict(task = task_1)
  Q0_0 <- sl_fit$predict(task = task_0)
  Q0_A <- sl_fit$predict(task = task)

  # Treatment probability Pr(A = 1)
  pi <- mean(data0[[treatment]] == 1)

  HA <- data0[[treatment]]/pi - (1 - data0[[treatment]])/(1 - pi)
  H1 <- dataA1[[treatment]]/pi
  H0 <- - (1 - dataA0[[treatment]])/(1 - pi)

  Q1 <- lm(data0[[primary]] ~ 0 + offset(Q0_A) + HA)
  epsilon <- coef(Q1)["HA"]

  Q1_1 <- Q0_1 + epsilon * H1
  Q1_0 <- Q0_0 + epsilon * H0
  Q1_A <- ifelse(data0[[treatment]], Q1_1, Q1_0)

  # Estimated treatment effect
  tmle <- mean(Q1_1 - Q1_0)

  # Influence function
  IC <- HA * (data0[[primary]] - Q1_A)
  sigma <- sqrt(var(IC)/nrow(data0))

  list(
    tmle_est = tmle,
    tmle_sigma = sigma,
    omega = sl_fit$coefficients[1],
    epsilon = epsilon,
    joint_sem = sl_fit$learner_fits$Lrnr_sem_NULL$fit_object
  )
}


#' Estimate an average treatment effect using model averaging (tmle3)
#' @inheritParams joint_sem
#' @examples
#' data(joint_example)
#' estimate_ma <- estimate_sem_ma_tmle(
#'   data0 = joint_example,
#'   endpoints = c("Y1", "Y2", "Y3"),
#'   categorical = c(),
#'   treatment = "A",
#'   primary = "Y1"
#' )
# estimate_sem_ma_tmle <- function(data0, endpoints, categorical = c(), treatment = "A", primary, ...) {
#
#   secondary <- endpoints[endpoints != primary]
#
#   # Set up tmle3
#   node_list <- list(
#     W = secondary,
#     A = treatment,
#     Y = primary
#   )
#   processed <- tmle3::process_missing(data0, node_list)
#   processed_data <- processed$data
#   node_list <- processed$node_list
#
#   ate_spec <- tmle3::tmle_ATE(
#     treatment_level = 1,
#     control_level = 0
#   )
#
#
#   # Chose base learners
#   lrnr_sem <- sl3::make_learner(Lrnr_sem, categorical = categorical)
#   lrnr_saturated <- sl3::make_learner(sl3::Lrnr_glm, covariates = treatment)
#
#   # Use mean learner for treatment model to get P(A) = (sum A)/n
#   lrnr_mean <- sl3::make_learner(sl3::Lrnr_mean)
#
#   # Super Learner for E(Y|A)
#   sl_Y <- sl3::Lrnr_sl$new(
#     learners = list(lrnr_sem, lrnr_saturated)
#   )
#
#   # Stack two mean learners for Pr(A) to prevent error with only one learner in
#   # the superlearner
#   sl_A <- sl3::Lrnr_sl$new(
#     learners = list(lrnr_mean, lrnr_mean)
#   )
#
#   learner_list <- list(A = sl_A, Y = sl_Y)
#
#   # TMLE ---
#   variable_type <- sl3:::variable_type()
#   tmle_fit <- tmle3::tmle3(ate_spec, processed_data, node_list, learner_list)
# }
