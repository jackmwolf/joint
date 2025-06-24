#' Verify input to [joint_sem()]
#'
#' @inheritParams joint_sem
#' @inheritParams ll_sem
#' @return TRUE if all input is valid
verify_input_joint_sem <- function(data0, endpoints, categorical, treatment) {

  # Verify input are characters
  if (!is.character(endpoints)) {
    stop("Input 'endpoints' must be a character vector.")
  }
  if (!is.character(categorical)) {
    stop("Input 'categorical' must be a character vector.")
  }
  if (!is.character(treatment)) {
    stop("Input 'treatment' must be a character vector.")
  }

  # Verify all endpoints are columns in data0
  if (any(!endpoints %in% colnames(data0))) {
    stop("Input in 'endpoints' does not correspond to a column in 'data0'.")
  }

  # Verify all categorical are listed in endpoints
  if (any(!categorical %in% endpoints)) {
    stop("Input in 'categorical' is not in 'endpoints'.")
  }

  # Verify all categorical are factors
  for (p in categorical) {
    if (!is.factor(data0[[p]])) stop(
      paste0(
        "All endpoints listed in 'categorical' must be factors. ",
        "Endpoint ",
        p,
        " is not a factor."
        )
    )
  }

  # Verify all not categorical are numeric
  for (p in setdiff(endpoints, categorical)) {
    if (!is.numeric(data0[[p]])) stop (
      paste0(
        "All endpoints not listed in 'categorical' must be numeric ",
        "Endpoint ",
        p,
        " is not numeric"
      )
    )
  }

  # Verify treatment is column in data0
  if (length(treatment) != 1) {
    stop("Input 'treatment' must be of length 1.")
  }
  if (!(treatment %in% colnames(data0))) {
    stop("Input in 'treatment' does not correspond to a column in 'data0'.")
  }

  # Verify >= 3 endpoints
  if (length(endpoints) < 3) {
    stop("Fewer than three endpoints input. joint_sem() requires 3 or more endpoints.")
  }
  return(TRUE)
}

#' Verify input to [joint_sl()]
#'
#' @inheritParams joint_sl
#' @return TRUE if all input is valid
verify_input_joint_sl <- function(data0, endpoints, treatment, primary, n_boot, ci_level) {

  # Check unique input for joint_sl
  if (!is.character(primary)) {
    stop("Input 'primary' must be a character vector.")
  }

  if (length(primary) != 1) {
    stop("Input 'primary' must be of length 1.")
  }

  if (!(primary %in% endpoints)) {
    stop("Input in 'primary' does not correpsond to a column in 'data0'.")
  }

  if ((n_boot < 0) | (n_boot != round(n_boot))) {
    stop("Inout 'n_boot' must be a non-negative interger")
  }

  if ((ci_level < 0) | (ci_level > 1)) {
    stop("Input 'ci_level' must be between 0 and 1")
  }

  # Check all other input used by joint_sem
  categorical <- endpoints[which(sapply(data0[, endpoints], is.factor))]
  verify_input_joint_sem(data0, endpoints, categorical, treatment)

  return(TRUE)
}
