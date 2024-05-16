# SUPER LEARNING HELPERS =======================================================
#' Wrapper function for \code{joint_sem} for Super Learning.
#' @param x Matrix of secondary endpoints and treatment indicator
#' @param y Vector of primary endpoint measurements
#' @param yname Character name for primary endpoint
#' @inheritParams joint_sem
#' @inherit joint_sem return
joint_sem_sl <- function(x, y, yname, categorical, ...) {

  data0 <- data.frame(y, x)
  names(data0)[1] <- yname
  endpoints <- colnames(data0)[colnames(data0) != "A"]

  fit_object <- joint_sem(
    data0 = data0, endpoints = endpoints, categorical = categorical)
  fit_object$A <- data0$A
  fit_object
}

#' Predict method for \code{joint_sem} for Super Learning.
#' @param object An object of class \code{joint_sem}
#' @param newdata An optional data frame with treatment indicators. If omitted
#'   the treatment indicators from the fitted model are used.
#' @return A vector of predicted outcomes given treatment indicators.
predict.joint_sem <- function(object, newdata) {
  yname <- object$endpoints[1]
  gamma <- object$estimate["gamma"]
  nu <- object$estimate[paste0("nu_", yname)]
  lambda <- object$estimate[paste0("lambda_", yname)]

  if (yname %in% object$categorical) {
    mu1 <- pnorm((nu + gamma * lambda)/sqrt(1 + lambda^2))
    mu0 <- pnorm((nu)/sqrt(1 + lambda^2))
  } else {
    mu1 <- nu + gamma * lambda
    mu0 <- nu
  }

  if (missing(newdata) || is.null(newdata)) {
    A <- object$A
  } else {
    A <- newdata$A
  }

  ifelse(A == 1, mu1, mu0)
}

##' \code{sl3} Learner for SEM \code{joint_sem}.
##'
##' @docType class
##' @importFrom R6 R6Class
##' @export
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
    .properties = c("continuous", "binomial"),

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
