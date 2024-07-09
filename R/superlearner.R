# SUPER LEARNING HELPERS =======================================================
#' Wrapper function for \code{joint_sem} for Super Learning.
#' @param x Matrix of secondary endpoints and treatment indicator
#' @param y Vector of primary endpoint measurements
#' @param yname Character name for primary endpoint
#' @param categorical Character vector of endpoint names to be treated as
#'   categorical (binary or ordinal) variables
#' @inheritParams joint_sem
#' @inherit joint_sem return
joint_sem_sl <- function(x, y, yname, categorical, ...) {

  data0 <- data.frame(y, x)

  # Undo temporary coding of factors as numeric to avoid sl3 creating
  # indicator variables.
  for (p in setdiff(categorical, yname)) {
    data0[[p]] <- factor(data0[[p]])
  }

  names(data0)[1] <- yname
  endpoints <- colnames(data0)[colnames(data0) != "A"]

  fit_object <- joint_sem(
    data0 = data0, endpoints = endpoints, ...)
  fit_object$A <- data0$A
  fit_object
}

#' Predict method for \code{joint_sem} for Super Learning.
#' @param object An object of class \code{joint_sem}
#' @param newdata An optional data frame with treatment indicators. If omitted
#'   the treatment indicators from the fitted model are used.
#' @importFrom stats pnorm
#' @importFrom sl3 pack_predictions
#' @return Predicted outcomes given treatment indicators.
predict.joint_sem <- function(object, newdata) {
  yname <- object$endpoints[1]
  gamma <- object$estimate["gamma"]
  nu <- object$estimate[paste0("nu_", yname)]
  lambda <- object$estimate[paste0("lambda_", yname)]

  if (missing(newdata) || is.null(newdata)) {
    A <- object$A
  } else {
    A <- newdata$A
  }

  # Ordinal endpoint: Matrix of predicted probabilities
  if (paste0("logdiffa_", yname, "_1") %in% names(object$estimate)) {

    # Conditional mean on link scale
    mu1 <- (nu + gamma * lambda)/sqrt(1 + lambda^2)
    mu0 <- (nu)/sqrt(1 + lambda^2)

    # Thresholds for integration
    logdiffs <- object$estimate[grep(paste0("logdiffa_", yname, "_"), names(object$estimate))]
    diffs <- exp(logdiffs)

    thresholds <- c(
      -Inf, 0,
      cumsum(diffs),
      Inf
    )
    names(thresholds) <- NULL

    # Predicted probabilities under A = 1
    probs_1 <- pnorm(thresholds, mean = mu1, sd = 1)
    preds_1 <- diff(probs_1)

    # Predicted probabilities under A = 0
    probs_0 <- pnorm(thresholds, mean = mu0, sd = 1)
    preds_0 <- diff(probs_0)

    # Return matrix of predicted probabilities
    out <- matrix(nrow = length(A), ncol = length(preds_1))
    out[A == 1, ] <- matrix(preds_1, nrow = sum(A == 1), ncol = length(preds_1), byrow = TRUE)
    out[A == 0, ] <- matrix(preds_0, nrow = sum(A == 0), ncol = length(preds_1), byrow = TRUE)

    # Return packed_predictions object for sl3 loss function
    out <- sl3::pack_predictions(out)

    return(out)

  } else {
    # Binary or numerical endpoint: Vector of expectations
    if (yname %in% object$categorical) {
      mu1 <- pnorm((nu + gamma * lambda)/sqrt(1 + lambda^2))
      mu0 <- pnorm((nu)/sqrt(1 + lambda^2))
    } else {
      mu1 <- nu + gamma * lambda
      mu0 <- nu
    }
    out <- ifelse(A == 1, mu1, mu0)

    return(out)
  }

}

##' \code{sl3} Learner for SEM \code{joint_sem}.
##'
##' @docType class
##' @importFrom R6 R6Class
##' @keywords data
##' @return Learner object with methods for training and prediction.
##'   See \code{\link{Lrnr_base}} for documentation on learners.
##' @format \code{\link{R6Class}} object.
##' @family Learners
##'
##' @section Parameters:
##' \describe{
##'   \item{\code{...}}{ Other parameters passed directly to \code{joint_sem}. See its documentation for details.
##'   }
##' }
##'
Lrnr_sem <- R6::R6Class(
  classname = "Lrnr_sem", inherit = sl3::Lrnr_base,
  portable = TRUE, class = TRUE,
  # Above, you should change Lrnr_template (in both the object name and the classname argument)
  # to a name that indicates what your learner does
  public = list(
    # you can define default parameter values here
    # if possible, your learner should define defaults for all required parameters
    # initialize = function(...) {
    #   super$initialize(params = sl3:::args_to_list(), ...)
    # }
  ),

  private = list(
    # list properties your learner supports here.
    # Use sl3_list_properties() for a list of options
    .properties = c("continuous", "binomial", "categorical"),

    # list any packages required for your learner here.
    .required_packages = c("joint"),

    # .train takes task data and returns a fit object that can be used to generate predictions
    .train = function(task) {
      # generate an argument list from the parameters that were
      # captured when your learner was initialized.
      # this allows users to pass arguments directly to your ml function
      args <- self$params

      # get outcome variable type
      # preferring learner$params$outcome_type first, then task$outcome_type
      outcome_type <- self$get_outcome_type(task)
      # should pass something on to your learner indicating outcome_type
      # e.g. family or objective

      # add task data to the argument list
      # what these arguments are called depends on the learner you are wrapping
      args$x <- as.matrix(task$X)
      args$y <- task$Y
      args$yname <- task$nodes$outcome

      # call a function that fits your algorithm
      # with the argument list you constructed

      fit_object <- sl3:::call_with_args(joint_sem_sl, args)

      # return the fit object, which will be stored
      # in a learner object and returned from the call
      # to learner$predict
      return(fit_object)
    },

    # .predict takes a task and returns predictions from that task
    .predict = function(task = NULL) {
      self$training_task
      self$training_outcome_type
      self$fit_object

      predictions <- predict(self$fit_object, task$X)
      return(predictions)
    }
  )
)
