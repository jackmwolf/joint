## FUNCTIONS FOR ESTIMATING SATURATED JOINT MODELS =============================

#' Estimate a saturated joint model
#'
#' `joint_saturated` fits a saturated joint model for all endpoints. Marginal
#' distributions are estimated through either linear or probit regression. The
#' appropriate variance and covariance terms are estimated through maximum
#' likelihood and the log likelihood is returned.
#'
#' @inheritParams joint_sem
#'
#' @inherit joint_sem return
#' @importFrom stats coef glm binomial lm residuals coefficients predict nlm
#' @importFrom VGAM vglm cumulative
joint_saturated <- function(data0, endpoints,
                            treatment = "A", ...) {
  t0 <- Sys.time()

  # List of categorical (factor) endpoints
  categorical <- endpoints[which(sapply(data0[, endpoints], is.factor))]

  verify_input_joint_sem(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment)

  n <- nrow(data0)
  P <- length(endpoints)
  q <- length(categorical)
  continuous <- endpoints[!(endpoints %in% categorical)]

  # Fit marginal models for each endpoint
  models <- vector(mode = "list", length = length(endpoints))
  names(models) <- endpoints

  # Store bounds for integration
  if (length(categorical) > 0) {
    lb <- matrix(nrow = n, ncol = length(endpoints))
    ub <- matrix(nrow = n, ncol = length(endpoints))
    colnames(lb) <- colnames(ub) <- endpoints
  }

  for (p in 1:length(endpoints)) {

    if (endpoints[p] %in% categorical) {

      if (length(unique(data0[, endpoints[p]])) > 2) {
        models[[p]] <- VGAM::vglm(
          as.numeric(data0[, endpoints[p]]) ~ data0[[treatment]],
          family = VGAM::cumulative(link="probitlink", parallel=TRUE)
          )
        thresholds <- c(
          -Inf,
          coef(models[[p]])[1:length(unique(data0[, endpoints[p]]))-1],
          Inf)
        lb[, p] <- thresholds[as.numeric(data0[, endpoints[p]])]
        ub[, p] <- thresholds[as.numeric(data0[, endpoints[p]]) + 1]

      } else {
        models[[p]] <- glm(
          data0[, endpoints[p]] ~ data0[[treatment]],
          family = binomial("probit"))
        thresholds <- c(-Inf, 0, Inf)
        lb[, p] <- thresholds[as.numeric(data0[, endpoints[p]])]
        ub[, p] <- thresholds[as.numeric(data0[, endpoints[p]]) + 1]

      }

    } else {
      models[[p]] <- lm(data0[, endpoints[p]] ~ data0[[treatment]])
    }
  }

  # Likelihood contributions ---
  ## Continuous endpoints ===
  if (P - q > 0) {
    # MLE of continuous by continuous endpoints
    E_continuous <- sapply(models[continuous], function(.x) .x$residuals)
    V_continuous <- 1/n * t(E_continuous) %*% E_continuous

    ll_continuous <- mvnfast::dmvn(
      X = E_continuous,
      mu = rep(0, ncol(E_continuous)),
      sigma = V_continuous,
      log = TRUE)

  } else {
    ll_continuous <- rep(0, n)
  }

  ## Categorical endpoints ===
  if (q > 0) {
    # Initialize categorical x categorical | continuous covariance
    E_categorical <- sapply(
      models[categorical],
      function(.x) {
        if (inherits(.x, "vglm")) {
          # residuals(.x, type = "response")[, 1]
          attributes(.x)$residuals[, 1]
        } else {
          residuals(.x, type = "deviance")
        }
      })

    E <- cbind(E_continuous, E_categorical)

    V_hat <- 1/n * t(E) %*% E
    scalar <- 1/sqrt(diag(V_hat))
    scalar[1:(P-q)] <- 1
    V_hat <- diag(scalar) %*% V_hat %*% diag(scalar)

    rownames(V_hat) <- colnames(V_hat) <- colnames(E)

    # Take required parameters of variance matrix for
    # Cov(X_continuous, X_categorical) and Cov(X_categorical, X_categorical)
    # by ignoring the first (P - q) * (P - q - 1)/2 elements of the upper
    # triangle
    n_exclude <- (P - q) * (P - q - 1)/2
    if (n_exclude == 0) {
      phi_init <- V_hat[upper.tri(V_hat)]
    } else {
      phi_init <- V_hat[upper.tri(V_hat)][-c(1:((P - q) * (P - q - 1)/2))]
    }

    mu <- sapply(
      models,
      function(.x) {
        if (inherits(.x, "vglm")) {
          data0[[treatment]] * coefficients(.x)[length(coefficients(.x))]
        } else if (inherits(.x, "glm")) {
          predict(.x, type = "link")
        } else if (inherits(.x, "lm")) {
          predict(.x)
        }
      })

    # Parameter estimation ---
    # Maximize -1 * log likelihood using nlm
    suppressWarnings({
      mle <- nlm(
        f = function(.x, ...) -1 * ll_cat_saturated(.x, ...),
        data0 = data0,
        endpoints = endpoints,
        categorical = categorical,
        treatment = treatment,
        continuous = continuous,
        lb = lb,
        ub = ub,
        mu = mu,
        Sigma = V_hat,
        p = phi_init,
        hessian = FALSE,
        stepmax = 2,
        ...
      )
    })
    ll_categorical <- -1 * mle$minimum
  } else {
    ll_categorical <- 0
  }

  t1 <- Sys.time()
  out <- list(
    # estimate = mle$estimate,
    ll = sum(ll_continuous) + ll_categorical,
    dim_phi = sum(sapply(models, function(.x) length(coef(.x)))) +
      P * (P - 1)/2 + (P - q),
    runtime = as.numeric(t1 - t0, units = "secs"),
    endpoints = endpoints,
    categorical = categorical
  )

  if (q > 0) {
    out$nlm_code <- mle$code
    out$nlm_iterations <- mle$iterations
  } else {
    out$nlm_iterations <- 0
  }

  return(out)

}

#' Calculate the log likelihood for categorical endpoints
#'
#' `ll_cat_saturated()` calculates the conditional log likelihood for all
#' categorical endpoints (conditional on all continuous endpoints) for input
#' marginal mean and variance matrices. As this function is used internally for
#' likelihood maximization, the primary input (`phi`) is a vector of covariance
#' terms for all continuous by categorical and categorical by categorical
#' covariances.
#'
#' @param phi Vector of covariance parameter estimates
#' @param continuous Vector of categorical endpoint names
#' @param lb,ub Matrix of lower and upper bounds for integration
#' @param mu Matrix of endpoint means (on the probit scale for categorical
#'   endpoints)
#' @param Sigma Matrix of endpoint covariances. Some entries are overwritten by
#'   `phi` before likelihood calculation.
#' @inheritParams joint_saturated
#' @inheritParams ll_sem
#' @importFrom MASS ginv
#' @importFrom mvtnorm pmvnorm
#' @return The log likelihood for the input parameters and data
ll_cat_saturated <- function(phi, data0, endpoints, categorical, treatment,
                             continuous, lb, ub, mu, Sigma) {

  P <- length(endpoints)
  q <- length(categorical)
  n <- nrow(data0)

  # V(Y) ---
  if (length(phi) == length(Sigma[upper.tri(Sigma)])) {
    Sigma[upper.tri(Sigma)] <- phi
  } else {
    Sigma[upper.tri(Sigma)][-c(1:((P - q) * (P - q - 1)/2))] <- phi
  }

  Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]

  colnames(Sigma) <- rownames(Sigma) <- c(continuous, categorical)

  ## Categorical endpoints ===
  if (P - q > 0) {

    inv_V22 <- MASS::ginv(Sigma[continuous, continuous, drop = FALSE])

    # conditional expectation
    E <- data0[, continuous, drop = FALSE] - mu[, continuous, drop = FALSE]
    c_e <- mu[, categorical, drop = FALSE] +
      t(Sigma[categorical, continuous, drop = FALSE] %*%
          inv_V22 %*%
          t(E)
      )

    # conditional variance
    c_v <- Sigma[categorical, categorical, drop = FALSE] -
      Sigma[categorical, continuous, drop = FALSE] %*%
      inv_V22 %*%
      Sigma[continuous, categorical, drop = FALSE]

  } else {
    c_e <- mu[, categorical, drop = FALSE]
    c_v <- Sigma[categorical, categorical, drop = FALSE]
  }

  ll_categorical <- vector(length = n)

  # Calculate probabilities for each observation by integrating over
  # q-dimensional hyper-rectangle bounded within lb[i, ], ub[i, ]
  for (i in 1:n) {
    ll_categorical[i] <- log(
      mvtnorm::pmvnorm(
        lower = lb[i, categorical],
        upper = ub[i, categorical],
        mean = c_e[i, categorical],
        sigma = c_v,
        keepAttr = FALSE))
  }

  return(sum(ll_categorical))
}
