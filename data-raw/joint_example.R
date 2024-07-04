## code to prepare `joint_example` dataset goes here

# Example data
set.seed(1)
n <- 200
A <- rep(c(0, 1), each = n/2)
nu0 <- rep(0, 4)
tau0 <- c(0.25, 0.5, -0.5, 0.25)
gamma0 <- 2

lambda0 <- tau0/gamma0
s0 <- 1 - (lambda0)^2
s0[3:4] <- 1

eta <- rnorm(n, mean = gamma0 * A)
EY <- rep(1, n) %*% t(nu0) + eta %*% t(lambda0)
e <- MASS::mvrnorm(n, mu = rep(0, length(nu0)), Sigma = diag(s0, nrow = 4, ncol = 4))
Y <- EY + e
Y3_cat <- ifelse(Y[, 3] > 0, 1, 0)
Y4_cat <- ifelse(Y[, 4] > qnorm(0.75), 2, ifelse(Y[, 4] > 0, 1, 0))

joint_example <- data.frame(A = A, Y1 = Y[, 1], Y2 = Y[, 2], Y3 = Y[, 3], Y4 = Y[, 4],
                    Y3_cat = factor(Y3_cat), Y4_cat = factor(Y4_cat))


usethis::use_data(joint_example, overwrite = TRUE)
