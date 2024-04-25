#' Verify input to joint_sem()
#'
#' @inheritParams joint_sem
#' @return TRUE if all input is valid
verify_input_joint_sem <- function(data0, endpoints, categorical, treatment) {

  # Verify all endpoints are columns in data0
  if (any(!endpoints %in% colnames(data0))) {
    stop("Input in 'endpoints' does not correspond to a column in 'data0'.")
  }

  # Verify all categorical are listed in endpoints
  if (any(!categorical %in% endpoints)) {
    stop("Input in 'categorical' is not in 'endpoints'.")
  }

  # Verify treatment is column in data0
  if (!(treatment %in% colnames(data0))) {
    stop("Input in 'treatment' does not correspond to a column in 'data0'.")
  }

  # Verify >= 3 endpoints
  if (length(endpoints) < 3) {
    stop("Fewer than three endpoints input. ate_sem() requires 3 or more endpoints.")
  }
  return(TRUE)
}
