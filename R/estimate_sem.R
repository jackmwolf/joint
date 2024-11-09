## FUNCTIONS FOR ESTIMATING JOINT MODELS WITH SEM CONSTRAINTS

#' Estimate a joint model using SEM constraints
#'
#' @param data0 Data frame with treatment assignment and endpoints
#' @param endpoints Character vector of endpoint names as in \code{data0}
#' @param treatment Character name of the column name corresponding to binary
#'   0-1 treatment indicator.
#' @param sandwich Logical indicator for whether to return the sandwich variance
#'   estimator in addition to the model-based estimator.
#' @param debug Logical indicator to print diagnostic information
#' @param phi_init Optional. Vector of initial parameter values for estimation.
#' @param ... additional arguments to [stats::nlm()].
#'
#' @importFrom numDeriv grad
#' @return An object of class [joint_sem]
#' @examples
#' data(joint_example)
#' fit <- joint_sem(
#'   joint_example, endpoints = c("Y1", "Y2", "Y3_cat"))
#' fit
#' estimate_effects(fit, risk_difference = "Y3_cat")
#' @export
joint_sem <- function(data0, endpoints, treatment = "A",
                      sandwich = FALSE, debug = FALSE,
                      phi_init = NULL, ...) {
  t0 <- Sys.time()

  # List of categorical (factor) endpoints
  categorical <- endpoints[which(sapply(data0[, endpoints], is.factor))]

  if (debug) cat("Verifying input\n")
  verify_input_joint_sem(
    data0 = data0, endpoints = endpoints, categorical = categorical,
    treatment = treatment)

  if (is.null(phi_init)) {
    if (debug) cat("Initalizing phi\n")
    phi_init <- initalize_phi(
      data0 = data0, endpoints = endpoints, categorical = categorical,
      treatment = treatment)
  } else {
    names(phi_init) <- names(
      initalize_phi(
        data0 = data0, endpoints = endpoints, categorical = categorical,
        treatment = treatment)
    )
  }
  if (debug) {
    cat("phi_init:\n")
    print(phi_init, digits = 2)
  }

  if (debug) cat("Maximizing log likelihood\n")
  # Parameter estimation ---

  # Maximize -1 * log likelihood using nlm
  suppressWarnings({
    mle <- nlm(
      f = function(.x, data0, endpoints, categorical, treatment, .names) {
        -1 * ll_sem(.x, data0, endpoints, categorical, treatment, .names)
      },
      data0 = data0,
      endpoints = endpoints,
      categorical = categorical,
      treatment = treatment,
      .names = names(phi_init),
      p = phi_init,
      hessian = TRUE,
      ...
    )
  })

  names(mle$estimate) <- names(phi_init)
  colnames(mle$hessian) <- rownames(mle$hessian) <- names(phi_init)

  if (debug) cat("Estimating model-based variance\n")
  # Variance estimation ---
  # E{H(theta)}^-1
  H_inv <- tryCatch(
    expr = {
      H_inv <- solve(mle$hessian)
      H_inv
    },
    error = function(cond) {
      message("Warning: joint_sem() could not invert hessian:")
      message(conditionMessage(cond))
      H_inv <- matrix(NA, ncol = ncol(mle$hessian), nrow = nrow(mle$hessian))
      rownames(H_inv) <- colnames(H_inv) <- rownames(mle$hessian)
      H_inv
    }
  )

  V_model <- H_inv
  V_out <- list(vcov = V_model)

  if (sandwich) {
    if (debug) cat("Estimating sandwich variance\n")
    # E{G(theta)^T(theta)}
    # Could improve computation by writing out closed form gradient/score function
    Gi <- sapply(
      split(data0, seq(nrow(data0))),
      function(.x) numDeriv::grad(
        ll_sem,
        x = mle$estimate,
        data0 = .x,
        endpoints = endpoints,
        categorical = categorical,
        treatment = treatment,
        .names = names(phi_init)
      )
    )

    GtG <- Gi %*% t(Gi)
    V_sand <- H_inv %*% GtG %*% H_inv
    dimnames(V_sand) <- dimnames(V_model)
    V_out$vcov_sandwich <- V_sand
  }

  # Output ---
  if (debug) cat("Preparing output\n")
  t1 <- Sys.time()
  out <- list(
    estimate = mle$estimate,
    ll = -1 * mle$minimum,
    dim_phi = length(mle$estimate),
    runtime = as.numeric(t1 - t0, units = "secs"),
    endpoints = endpoints,
    categorical = categorical
  )

  out <- append(out, V_out)

  nlm_diagnostic <- list(
    nlm_code = mle$code,
    nlm_iterations = mle$iterations
  )
  out <- append(out, nlm_diagnostic)

  class(out) <- "joint_sem"
  return(out)
}

#' Create initial vector of parameter values for estimating a joint model
#'
#' @inheritParams joint_sem
#' @param categorical Character vector of endpoint names to be treated as
#'   categorical (binary or ordinal) variables
#' @param gamma Initial value of the gamma parameter, defaults to 2
#' @importFrom stats qnorm var
#' @return A named vector of initial parameter values for likelihood
#'   maximization.
initalize_phi <- function(data0, endpoints, categorical = c(), treatment, gamma = 2) {
  # Initialize parameters ---
  P <- length(endpoints)
  q <- length(categorical)

  if (q > 0) {
    # Identify appropriate # levels for categorical endpoints
    n_levels <- vector(mode = "integer", length = length(categorical))
    for (i in seq_along(categorical)) {
      n_levels[i] <- length(levels(data0[[categorical[i]]]))
    }
    names(n_levels) <- categorical

    # Number of required threshold parameters for partitioning real line for
    # the probit model (-infty < 0 < a1 < a2 < ... < a_nthresholds < infty)
    n_thresholds <- n_levels - 2
  }


  nu <- vector(mode = "double", length = P)
  lambda <- vector(mode = "double", length = P)
  logtheta <- vector(mode = "double", length = 0)
  a <- vector(mode = "double", length = 0)
  logdiffa <- vector(mode = "double", length = 0)


  for (p in 1:length(endpoints)) {
    if (endpoints[p] %in% categorical) {

      nu[p] <- qnorm(
        1 - mean(
          as.numeric(data0[data0[[treatment]] == 0, endpoints[p]]) == 1
          )
        )
      lambda[p] <- 1/gamma * (
        1 - mean(
          as.numeric(data0[data0[[treatment]] == 1, endpoints[p]]) == 1
          ) - nu[p])

      if (n_thresholds[endpoints[p]] > 0) {

        cum_props_p <- cumsum(prop.table(table(data0[data0[[treatment]] == 0, endpoints[p]])))
        a_p <- qnorm(cum_props_p)[2:(length(cum_props_p) - 1)] - qnorm(cum_props_p)[1]
        logdiffa_p <- log(diff(c(0, a_p)))
        names(a_p) <- paste0("a_", endpoints[p], "_", 1:length(a_p))
        names(logdiffa_p) <- paste0("logdiffa_", endpoints[p], "_", 1:length(a_p))
        a <- c(a, a_p)
        logdiffa <- c(logdiffa, logdiffa_p)

      }

    } else {

      nu[p] <- mean(data0[data0[[treatment]] == 0, endpoints[p]])
      lambda[p] <- 1/gamma * (mean(data0[data0[[treatment]] == 1, endpoints[p]]) - nu[p])
      logtheta_p <- log(max(var(data0[data0[[treatment]] == 0, endpoints[p]]) - lambda[p]^2, 1e-5))
      names(logtheta_p) <- paste0("logtheta_", endpoints[p])
      logtheta <- c(logtheta, logtheta_p)

    }

    names(nu)[p] <- paste0("nu_", endpoints[p])
    names(lambda)[p] <- paste0("lambda_", endpoints[p])
  }

  phi_init <- c(
    gamma,
    nu,
    lambda,
    logtheta,
    logdiffa
  )

  names(phi_init)[1] <- "gamma"
  return(phi_init)
}

#' Calculate the log likelihood for a joint model using SEM constraints
#'
#' @param phi Named vector of parameter values
#' @param .names Optimal vector of names for phi
#' @inheritParams joint_sem
#' @param categorical Character vector of endpoint names to be treated as
#'   categorical (binary or ordinal) variables
#' @importFrom mvnfast dmvn
#' @return The log likelihood for the input parameters and data
ll_sem <- function(phi, data0, endpoints, categorical, treatment, .names = NULL) {

  if (!is.null(.names)) names(phi) <- .names

  P <- length(endpoints)
  q <- length(categorical)

  continuous <- endpoints[!(endpoints %in% categorical)]

  gamma <- phi["gamma"]
  nu <- phi[paste0("nu_", endpoints)]
  lambda <- phi[paste0("lambda_", endpoints)]
  theta <- vector(mode = "double", length = P)

  for (p in 1:length(endpoints)) {
    if (endpoints[p] %in% categorical) {
      theta[p] <- 1
    } else {
      theta[p] <- exp(phi[paste0("logtheta_", endpoints[p])])
    }
  }
  names(theta) <- paste0("theta_", endpoints)

  # if (any(theta < 0)) return(Inf)

  n <- nrow(data0)

  # E(Y, Y*) ---
  mu <- matrix(1, n) %*% matrix(nu, nrow = 1) +
    as.matrix(data0[[treatment]], ncol = 1) %*% matrix(gamma * lambda, nrow = 1)
  colnames(mu) <- endpoints

  # V(Y) ---
  Sigma <- lambda %*% t(lambda) + diag(theta, nrow = P, ncol = P)
  colnames(Sigma) <- rownames(Sigma) <- endpoints
  if (any(eigen(Sigma)$values <= 0)) return(Inf)

  # Likelihood contributions ---
  ## Continuous endpoints ===
  if (P - q > 0) {
    E <- data0[, continuous, drop = FALSE] - mu[, continuous, drop = FALSE]
    V <- Sigma[continuous, continuous, drop = FALSE]

    ll_continuous <- tryCatch(
      expr = {
        mvnfast::dmvn(X = as.matrix(E), mu = rep(0, ncol(E)), sigma = V,
                      log = TRUE)
      },
      error = function(e) {
        # message("Error: ", e$message)
        NA
      }
    )

  } else {
    ll_continuous <- rep(0, n)
  }

  ## Categorical endpoints ===
  if (q > 0) {

    if (P - q > 0) {
      # Conditional mean and variance
      # conditional expectation
      # print(Sigma[categorical, categorical, drop = FALSE])
      c_e <- mu[, categorical, drop = FALSE] +
        t(Sigma[categorical, continuous, drop = FALSE] %*%
            ginv(Sigma[continuous, continuous, drop = FALSE]) %*%
            t(E)
        )

      # conditional variance
      c_v <- Sigma[categorical, categorical, drop = FALSE] -
        Sigma[categorical, continuous, drop = FALSE] %*%
        ginv(V[continuous, continuous, drop = FALSE]) %*%
        Sigma[continuous, categorical, drop = FALSE]

    } else {
      c_e <- mu[, categorical, drop = FALSE]
      c_v <- Sigma[categorical, categorical, drop = FALSE]
    }

    # Set up bounds for integration
    lb <- matrix(nrow = n, ncol = length(categorical))
    ub <- matrix(nrow = n, ncol = length(categorical))

    for (p in 1:length(categorical)) {
      logdiffs <- phi[grep(paste0("logdiffa_", categorical[p], "_"), names(phi))]
      diffs <- exp(logdiffs)

      thresholds <- c(
        -Inf, 0,
        cumsum(diffs),
        Inf
      )
      lb[, p] <- thresholds[as.numeric(data0[, categorical[p]])]
      ub[, p] <- thresholds[as.numeric(data0[, categorical[p]]) + 1]
    }

    ll_categorical <- vector(length = n)

    # Calculate probabilities for each observation by integrating over
    # q-dimensional hyper-rectangle bounded within lb[i, ], ub[i, ]
    # SLOW AND REQUIRES INVERTING SIGMA EVERY TIME
    for (i in 1:n) {
      ll_categorical[i] <- log(
        mvtnorm::pmvnorm(
          lower = lb[i, ], upper = ub[i, ], mean = c_e[i, ], sigma = c_v,
          keepAttr = FALSE))
    }
  } else {
    ll_categorical <- rep(0, n)
  }

  return(sum(ll_continuous + ll_categorical))
}


#' Estimate treatment effects from a joint model
#'
#' @param model An object of class [joint_sem]
#' @param risk_difference Names of endpoints for which to calculate the risk
#'   difference.
#' @param sandwich Logical indicator to use the sandwich variance estimator
#' @inherit joint_sem examples
#' @importFrom stats pnorm dnorm
#' @return A list of effect estimates and variances
#' @export
estimate_effects <- function(model, risk_difference = c(), sandwich = FALSE) {

  estimate <- model$estimate

  if (sandwich) {
    if (!is.null(model$vcov_sandwich)) {
      vcov <- model$vcov_sandwich
    } else {
      warning("Input model does not have sandwich variance matrix. Using model-based variance.")
      vcov <- model$vcov
    }
  } else {
    vcov <- model$vcov
  }

  endpoints <- model$endpoints

  params <- c("gamma", paste0("lambda_", endpoints))
  if (length(risk_difference) > 0) {
    params <- c(params, paste0("nu_", risk_difference))
  }

  beta <- estimate[params]
  V <- vcov[params, params]

  # Vector for storing point estimates
  ests <- vector(mode = "double", length = length(endpoints) + length(risk_difference) * 2)
  if (length(risk_difference) > 0) {
    names(ests) <- c(endpoints, paste0(risk_difference, "_RD"), paste0(risk_difference, "_SD"))
  } else {
    names(ests) <- endpoints
  }


  # Matrix of gradient for transformation beta --> ests
  j <- matrix(0, nrow = length(ests), ncol = length(beta))
  rownames(j) <- names(ests)
  colnames(j) <- names(beta)

  # Effects on link scale
  for (y in endpoints) {
    ests[y] <- beta["gamma"] * beta[paste0("lambda_", y)]
    j[y, "gamma"] <- beta[paste0("lambda_", y)]
    j[y, paste0("lambda_", y)] <- beta["gamma"]
  }

  # Effects on risk difference scale and probit coefficients (with var = 1)
  for (y in risk_difference) {
    gamma <- beta["gamma"]
    nu <- beta[paste0("nu_", y)]
    lambda <- estimate[paste0("lambda_", y)]

    ests[paste0(y, "_RD")] <- unname(
      pnorm((nu + gamma * lambda)/sqrt(1 + lambda^2)) -
        pnorm((nu)/sqrt(1 + lambda^2))
    )

    d1 <- dnorm((nu + gamma * lambda)/sqrt(1 + lambda^2))
    d0 <- dnorm((nu)/sqrt(1 + lambda^2))

    j[paste0(y, "_RD"), "gamma"] <- d1 * lambda/sqrt(1 + lambda^2)
    j[paste0(y, "_RD"), paste0("nu_", y)] <- (d1 - d0) * 1/sqrt(1 + lambda^2)
    j[paste0(y, "_RD"), paste0("lambda_", y)] <-
      (d1 * (gamma - nu * lambda) - d0 * (-nu * lambda))/
      (1 + lambda^2)^(3/2)

    ests[paste0(y, "_SD")] <- unname((gamma * lambda)/sqrt(1 + lambda^2))
    j[paste0(y, "_SD"), "gamma"] <- lambda/sqrt(1 + lambda^2)
    j[paste0(y, "_SD"), paste0("lambda_", y)] <- gamma/sqrt(1 + lambda^2)

  }

  v_est <- j %*% V %*% t(j)

  out <- list(
    estimate = ests,
    vcov = v_est
  )
  return(out)
}

estimate_marginal <- function(model) {

  estimate <- model$estimate
  endpoints <- model$endpoints

  gamma <- estimate["gamma"]
  nu <- estimate[paste0("nu_", endpoints)]
  lambda <- estimate[paste0("lambda_", endpoints)]

  theta <- rep(1, length = length(endpoints))
  names(theta) <- paste0("theta_", endpoints)

  for (y in endpoints) {
    if (paste0("logtheta_", y) %in% names(estimate)) {
      theta[paste0("theta_", y)] <- exp(estimate[paste0("logtheta_", y)])
    }
  }

  mu0 <- nu
  mu1 <- nu + gamma * lambda

  names(mu0) <- names(mu1) <- endpoints

  V <- lambda %*% t(lambda) + diag(theta)

  out <- list(
    mu0 = mu0,
    mu1 = mu1,
    V = V
  )
  return(out)
}

